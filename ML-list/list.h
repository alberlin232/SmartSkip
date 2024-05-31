/*
 * File:
 *   skiplist.h
 * Author(s):
 *   Vincent Gramoli <vincent.gramoli@epfl.ch>
 * Description:
 *   Stress test of the skip list implementation.
 *
 * Copyright (c) 2009-2010.
 *
 * skiplist.h is part of Synchrobench
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

#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <pthread.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>

#include <atomic_ops.h>

#include "radix_spline.h"

extern volatile AO_t stop;
extern unsigned int global_seed;
#ifdef TLS
extern __thread unsigned int *rng_seed;
#else /* ! TLS */
extern pthread_key_t rng_seed_key;
#endif /* ! TLS */
extern unsigned int levelmax;

#define TRANSACTIONAL                   d->unit_tx

typedef intptr_t val_t;
typedef intptr_t level_t;
#define VAL_MIN                         INT_MIN
#define VAL_MAX                         INT_MAX


typedef struct sl_node {
  val_t val;
  bool deleted;
  struct sl_node *next;
} sl_node_t;

typedef struct sl_intset {
  sl_node_t *head;
} sl_intset_t;

typedef struct shift_node {
  int count;
  int delta;
  sl_node_t* node;
}shift_node_t;

sl_node_t *sl_new_node(val_t val, sl_node_t *next, int transactional);
void sl_delete_node(sl_node_t *n);

sl_intset_t *sl_set_new();
void sl_set_delete(sl_intset_t *set);
unsigned long sl_set_size(sl_intset_t *set);

shift_node_t *new_shift_table(int table_size);
void populate_shift_table(sl_intset_t *set, shift_node_t *shift_table, RadixSpline<val_t> *spline, int table_size);

int seq_add(sl_intset_t *set, val_t val);
int sl_contains(sl_intset_t *set, RadixSpline<val_t> *spline, shift_node_t *shift_table, int table_size, val_t val, unsigned long *iterations);
int sl_add(sl_intset_t *set, RadixSpline<val_t> *spline, shift_node_t *shift_table, int table_size, val_t val, unsigned long *iterations);
int sl_remove(sl_intset_t *set, RadixSpline<val_t> *spline, shift_node_t *shift_table, int table_size, val_t val, unsigned long *iterations);
