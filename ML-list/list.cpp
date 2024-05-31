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

#include "list.h"	


/*
 * Returns a random level for inserting a new node, results are hardwired to p=0.5, min=1, max=32.
 *
 * "Xorshift generators are extremely fast non-cryptographically-secure random number generators on
 * modern architectures."
 *
 * Marsaglia, George, (July 2003), "Xorshift RNGs", Journal of Statistical Software 8 (14)
 */

/* 
 * Create a new node with its next field. 
 * If next=NULL, then this create a tail node. 
 */
sl_node_t *sl_new_node(val_t val, sl_node_t *next, int transactional)
{
  sl_node_t *node;

  node = (sl_node_t *)malloc(sizeof(sl_node_t));
  if (node == NULL) {
	perror("malloc");
	exit(1);
  }

  node->val = val;
  node->next = next;
  node->deleted = false;

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
  max = sl_new_node(VAL_MAX, NULL, 0);
  min = sl_new_node(VAL_MIN, max, 0);
  set->head = min;

  return set;
}

void sl_set_delete(sl_intset_t *set)
{
  sl_node_t *node, *next;

  node = set->head;
  while (node != NULL) {
    next = node->next;
    sl_delete_node(node);
    node = next;
  }
  free(set);
}

unsigned long sl_set_size(sl_intset_t *set)
{
  unsigned long size = 0;
  sl_node_t *node;

  node = set->head->next;
  while (node->next != NULL) {
    size++;
    node = node->next;
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
  while (node->next->next != NULL) {
    node = node->next;
    int k = spline->GetEstimatedPosition(node->val) * (table_size-1);
    int delta = j - k;
    if (delta <= shift_table[k].delta) {
      shift_table[k].delta = delta;
      shift_table[k].node = node;
    }
    shift_table[k].count++;
    j++;
  }


  for (j = table_size-1; j >= 0; j--) {
    if (shift_table[j].count == 0) {
      shift_table[j].count = shift_table[j+1].count;
      shift_table[j].delta = shift_table[j+1].delta+1;
      shift_table[j].node = shift_table[j+1].node;
    }
  }
}

int sl_contains(sl_intset_t *set, RadixSpline<val_t> *spline, shift_node_t *shift_table, int table_size, val_t val, unsigned long *iterations)
{
	int result = 0;
	sl_node_t *node;
	
	int k = spline->GetEstimatedPosition(val) * (table_size-1);
  node = shift_table[k].node;

  while(node->val > val && k-- > 0) {
    node = shift_table[k].node;
  }
  
  while (node->val < val) {
    (*iterations)++;
    node = node->next;
  }
  if (node->val == val && !node->deleted) {
    result = 1;
  }

	return result;
}


int sl_add(sl_intset_t *set, RadixSpline<val_t> *spline, shift_node_t *shift_table, int table_size, val_t val, unsigned long *iterations) {
  sl_node_t *pred, *curr = NULL;

  int k = spline->GetEstimatedPosition(val) * (table_size-1);
  pred = shift_table[k].node;

  while(pred->val > val && k-- > 0) {
    pred = shift_table[k].node;
  }

  pred = set->head;
  curr = pred->next;
  while (curr->val < val) {
    pred = curr;
    curr = curr->next;
    (*iterations)++;
  }

  if (curr->val == val ) {
    if (curr->deleted) {
      curr->deleted = false;
      return 1;
    }
    return 0;
  }

  sl_node *newnode = sl_new_node(val, curr, 0);
	pred->next = newnode;

	return 1;
}

int sl_remove(sl_intset_t *set, RadixSpline<val_t> *spline, shift_node_t *shift_table, int table_size, val_t val, unsigned long *iterations) {
  sl_node_t *pred, *curr = NULL;

  int k = spline->GetEstimatedPosition(val) * (table_size-1);
  pred = shift_table[k].node;

  while(pred->val > val && k-- > 0) {
    pred = shift_table[k].node;
  }

  curr = pred->next;
  while (curr->val < val) {
    pred = curr;
    curr = curr->next;
    (*iterations)++;
  }

  if (curr->val != val || curr->deleted) {
    return 0;
  }
  curr->deleted = true;

	return 1;
}

int seq_add(sl_intset_t *set, val_t val) {
  sl_node_t *pred, *curr = NULL;

  pred = set->head;
  curr = pred->next;
  while (curr->val < val) {
    pred = curr;
    curr = curr->next;
  }

  if (curr->val == val) {
    return 0;
  }
  sl_node *newnode = sl_new_node(val, curr, 0);
	pred->next = newnode;

	return 1;
}