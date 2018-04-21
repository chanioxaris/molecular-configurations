#ifndef _FUNCTIONS_
#define _FUNCTIONS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <cblas.h>
#include <lapacke.h>

#include "structs.c"

#define PATH_LENGTH 32
#define LINE_LENGTH 26000

#define THRESHOLD 0.05


// Declare functions of functions.c file
int is_centroid(int *, int, int);

double rand_gaussian();

int centroids_transposition(conformation **, int *, int *, int, double **);

double **create_distances_array(conformation **, int, int, int);

#endif
