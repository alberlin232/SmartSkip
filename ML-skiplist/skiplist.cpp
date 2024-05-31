/*
 * File:
 *   skiplist.c
 * Author(s):
 *   Vincent Gramoli <vincent.gramoli@epfl.ch>
 * Description:
 *   Skip list definition 
 *
 * Copyright (c) 2009-2010.
 *
 * skiplist.c is part of Synchrobench
 * 
 * Synchrobench is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, version 2
 * of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include "skiplist.h"	

unsigned int levelmax = MAXLEVEL;

/*
 * Returns a random level for inserting a new node, results are hardwired to p=0.5, min=1, max=32.
 *
 * "Xorshift generators are extremely fast non-cryptographically-secure random number generators on
 * modern architectures."
 *
 * Marsaglia, George, (July 2003), "Xorshift RNGs", Journal of Statistical Software 8 (14)
 */
int get_rand_level() {
	static uint32_t y = 2463534242UL;
	y^=(y<<13);
	y^=(y>>17);
	y^=(y<<5);
	uint32_t temp = y;
	uint32_t level = 1;
	while (((temp >>= 1) & 1) != 0) {
		++level;
	}
	/* 1 <= level <= levelmax */
	if (level > levelmax) {
		return (int)levelmax;
	} else {
		return (int)level;
	}
}

int floor_log_2(unsigned int n) {
  int pos = 0;
  if (n >= 1<<16) { n >>= 16; pos += 16; }
  if (n >= 1<< 8) { n >>=  8; pos +=  8; }
  if (n >= 1<< 4) { n >>=  4; pos +=  4; }
  if (n >= 1<< 2) { n >>=  2; pos +=  2; }
  if (n >= 1<< 1) {           pos +=  1; }
  return ((n == 0) ? (-1) : pos);
}

/* 
 * Create a new node without setting its next fields. 
 */
sl_node_t *sl_new_simple_node(val_t val, int toplevel, int transactional)
{
  sl_node_t *node;
  node = (sl_node_t *)malloc(sizeof(sl_node_t) + toplevel * sizeof(sl_node_t *));
  if (node == NULL) {
    perror("malloc");
    exit(1);
    }

  node->val = val;
  node->toplevel = toplevel;
  node->deleted = 0;

  return node;
}

/* 
 * Create a new node with its next field. 
 * If next=NULL, then this create a tail node. 
 */
sl_node_t *sl_new_node(val_t val, sl_node_t *next, int toplevel, int transactional)
{
  sl_node_t *node;
  int i;

  node = sl_new_simple_node(val, toplevel, transactional);

  for (i = 0; i < levelmax; i++)
    node->next[i] = next;
	
  return node;
}

void sl_delete_node(sl_node_t *n)
{
  free(n);
}

sl_intset_t *sl_set_new()
{
  sl_intset_t *set;
  sl_node_t *min, *max;
	
  if ((set = (sl_intset_t *)malloc(sizeof(sl_intset_t))) == NULL) {
    perror("malloc");
    exit(1);
  }
  max = sl_new_node(VAL_MAX, NULL, levelmax, 0);
  min = sl_new_node(VAL_MIN, max, levelmax, 0);
  set->head = min;
  return set;
}

void sl_set_delete(sl_intset_t *set)
{
  sl_node_t *node, *next;

  node = set->head;
  while (node != NULL) {
    next = node->next[0];
    sl_delete_node(node);
    node = next;
  }
  free(set);
}

unsigned long sl_set_size(sl_intset_t *set)
{
  unsigned long size = 0;
  sl_node_t *node;

  node = set->head->next[0];
  while (node->next[0] != NULL) {
    if (!node->deleted)
      size++;
    node = node->next[0];
  }

  return size;
}

shift_node_t *new_shift_table(int table_size) {
  shift_node_t *shift_table = (shift_node_t *) malloc((table_size) * sizeof(shift_node_t));
  for (int i = 0; i < table_size; i++) {
    shift_table[i].count = 0;
    shift_table[i].delta = INT_MAX;
    shift_table[i].node = NULL;
  }
  return shift_table;
}

void populate_shift_table(sl_intset_t *set, shift_node_t *shift_table, RadixSpline<val_t> *spline, int table_size) {
  sl_node_t *node = set->head;
  int j = 0;

  //put head at index 0
  shift_table[0].node = node;
  shift_table[0].count = 1;
  shift_table[0].delta = 0;

  while (node->next[0]->next[0] != NULL) {
    node = node->next[0];
    // if (node->toplevel != levelmax)
    //   continue;
    int k = spline->GetEstimatedPosition(node->val) * (table_size-1);
    int delta = j - k;
    if (delta <= shift_table[k].delta) {
      shift_table[k].delta = delta;
      shift_table[k].node = node;
    }
    shift_table[k].count++;
    j++;
  }

  //put tail at index table_size-1
  shift_table[table_size-1].node = node->next[0];
  shift_table[table_size-1].count = 1;
  shift_table[table_size-1].delta = 0;


  for (j = table_size-1; j >= 0; j--) {
    if (shift_table[j].count == 0) {
      shift_table[j].count = shift_table[j+1].count;
      shift_table[j].delta = shift_table[j+1].delta+1;
      shift_table[j].node = shift_table[j+1].node;
    }
  }
}

inline int is_marked(uintptr_t i) {
  return (int)(i & (uintptr_t)0x01);
}

inline uintptr_t unset_mark(uintptr_t i) {
  return (i & ~(uintptr_t)0x01);
}

inline uintptr_t set_mark(uintptr_t i) {
  return (i | (uintptr_t)0x01);
}

inline void fraser_search(sl_intset_t *set, 
                          val_t val, 
                          sl_node_t **left_list, 
                          sl_node_t **right_list,
                          RadixSpline<val_t> *spline, 
                          shift_node_t *shift_table, 
                          int table_size, 
                          unsigned long *iterations) {
  int i;
  sl_node_t *left, *left_next, *right, *right_next;

retry:

  int k = spline->GetEstimatedPosition(val) * (table_size-1);
  left = (sl_node_t *) unset_mark((long) shift_table[k].node);
  (*iterations)++;
  while ((left->val > val && k-- > 0)) {
    left = (sl_node_t *) unset_mark((long) shift_table[k].node);
    (*iterations)++;
  }
  // left = (sl_node_t *) unset_mark((long) left->next);

  for (i = left->toplevel - 1; i >= 0; i--) {
    left_next = left->next[i];
    if (is_marked((uintptr_t)left_next))
      goto retry;
    /* Find unmarked node pair at this level */
    for (right = left_next; ; right = right_next) {
      /* Skip a sequence of marked nodes */
      while(1) {
        right_next = right->next[i];
        if (!is_marked((uintptr_t)right_next))
          break;
        right = (sl_node_t*)unset_mark((uintptr_t)right_next);
      }
      if (right->val >= val)
        break;
      left = right; 
      left_next = right_next;
    }
    /* Ensure left and right nodes are adjacent */
    if ((left_next != right) && 
        (!ATOMIC_CAS_MB(&left->next[i], left_next, right)))
      goto retry;
    if (left_list != NULL)
      left_list[i] = left;
    if (right_list != NULL)	
      right_list[i] = right;
  }
}

inline void mark_node_ptrs(sl_node_t *n) {
  int i;
  sl_node_t *n_next;
	
  for (i=n->toplevel-1; i>=0; i--) {
    do {
      n_next = n->next[i];
      if (is_marked((uintptr_t)n_next))
        break;
    } while (!ATOMIC_CAS_MB(&n->next[i], n_next, set_mark((uintptr_t)n_next)));
  }
}

int sl_contains(sl_intset_t *set, 
                RadixSpline<val_t> *spline, 
                shift_node_t *shift_table, 
                int table_size, 
                val_t val, 
                unsigned long *iterations)
{
  sl_node_t **succs;
  int result;

  succs = (sl_node_t **)malloc(levelmax * sizeof(sl_node_t *));
  fraser_search(set, val, NULL, succs, spline, shift_table, table_size, iterations);
  result = (succs[0]->val == val && !succs[0]->deleted);
  free(succs);
  return result;
}


int sl_add(sl_intset_t *set, 
           RadixSpline<val_t> *spline, 
           shift_node_t *shift_table, 
           int table_size, val_t v, 
           unsigned long *iterations) 
{
  sl_node_t *new_n, *new_next, *pred, *succ, **succs, **preds;
  int i;
  int result;

  new_n = sl_new_simple_node(v, get_rand_level(), 6);
  preds = (sl_node_t **)malloc(levelmax * sizeof(sl_node_t *));
  succs = (sl_node_t **)malloc(levelmax * sizeof(sl_node_t *));
retry: 	
  fraser_search(set, v, preds, succs, spline, shift_table, table_size, iterations);
  /* Update the value field of an existing node */
  if (succs[0]->val == v) {
    /* Value already in list */
    if (succs[0]->deleted) {
      /* Value is deleted: remove it and retry */
      mark_node_ptrs(succs[0]);
      goto retry;
    }
    result = 0;
    sl_delete_node(new_n);
    goto end;
  }
  for (i = 0; i < new_n->toplevel; i++)
    new_n->next[i] = succs[i];
  /* Node is visible once inserted at lowest level */
  if (!ATOMIC_CAS_MB(&preds[0]->next[0], succs[0], new_n)) 
    goto retry;
  for (i = 1; i < new_n->toplevel; i++) {
    while (1) {
      pred = preds[i];
      succ = succs[i];
      /* Update the forward pointer if it is stale */
      new_next = new_n->next[i];
      if ((new_next != succ) && 
          (!ATOMIC_CAS_MB(&new_n->next[i], unset_mark((uintptr_t)new_next), succ)))
        break; /* Give up if pointer is marked */
      /* Check for old reference to a k node */
      if (succ->val == v)
        succ = (sl_node_t *)unset_mark((uintptr_t)succ->next);
      /* We retry the search if the CAS fails */
      if (ATOMIC_CAS_MB(&pred->next[i], succ, new_n))
        break;
      fraser_search(set, v, preds, succs, spline, shift_table, table_size, iterations);
    }
  }
  result = 1;
end:
  free(preds);
  free(succs);

  return result;
}

int sl_remove(sl_intset_t *set, 
              RadixSpline<val_t> *spline, 
              shift_node_t *shift_table, 
              int table_size, 
              val_t val, 
              unsigned long *iterations)
{
 sl_node_t **succs;
  int result;

  succs = (sl_node_t **)malloc(levelmax * sizeof(sl_node_t *));
  fraser_search(set, val, NULL, succs, spline, shift_table, table_size, iterations);
  result = (succs[0]->val == val);
  if (result == 0)
    goto end;
  /* 1. Node is logically deleted when the deleted field is not 0 */
  if (succs[0]->deleted) {
    result = 0;
    goto end;
  }
  ATOMIC_FETCH_AND_INC_FULL(&succs[0]->deleted);
  /* 2. Mark forward pointers, then search will remove the node */
  mark_node_ptrs(succs[0]);
  fraser_search(set, val, NULL, NULL, spline, shift_table, table_size, iterations);    
end:
  free(succs);

  return result;
}

int seq_add(sl_intset_t *set, val_t val) {
	int i, l, result;
	sl_node_t *node, *next;
	sl_node_t *preds[MAXLEVEL], *succs[MAXLEVEL];
	
	node = set->head;
	for (i = node->toplevel-1; i >= 0; i--) {
		next = node->next[i];
		while (next->val < val) {
			node = next;
			next = node->next[i];
		}
		preds[i] = node;
		succs[i] = node->next[i];
	}
	node = node->next[0];
	if ((result = (node->val != val)) == 1) {
		l = get_rand_level();
		node = sl_new_simple_node(val, l, 0);
		for (i = 0; i < l; i++) {
			node->next[i] = succs[i];
			preds[i]->next[i] = node;
		}
	}
	return result;
}