#include "functions.h"
#include "metric_functions.h"

// Function that calculates the Frechet distance between two curves (using a 2D array)
double frechet_distance(double **conf1, double **conf2, int N) {
  	int i, j;
	double distance;
	double **array;

  	array = (double**) malloc(N * sizeof(double*));

	if (array == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

    for (i = 0 ; i < N ; i++) {
    	array[i] = (double*) malloc(N * sizeof(double));

		if (array[i] == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}

    array[0][0] = euclidean(conf1[0], conf2[0]);

    for (i = 1 ; i < N ; i++)
    	array[i][0] = max_2(array[i-1][0], euclidean(conf1[i], conf2[0]));

    for (j = 1 ; j < N ; j++)
    	array[0][j] = max_2(array[0][j-1], euclidean(conf1[0], conf2[j]));

    for (i = 1 ; i < N ; i++)
      	for (j = 1 ; j < N ; j++)
        	array[i][j] = max_2(min_3(array[i-1][j], array[i][j-1], array[i-1][j-1]), euclidean(conf1[i], conf2[j]));

	distance = array[N-1][N-1];

    for (i = 0 ; i < N ; i++)
    	free(array[i]);
    free(array);

	return distance;
}


// Function that finds out the max value between two numbers
double max_2(double a, double b) {
	if (a >= b)
    	return a;

  	return b;
}


// Function that finds out the min value out of three numbers
double min_3(double a, double b, double c) {
	if (a <= b && a <= c)
      	return a;
  	else if (b <= a && b <= c)
      	return b;
    else
      	return c;
}


// Function that calculates the Euclidean distance between two points
double euclidean(double *coord1, double *coord2) {
	int i;
  	double dist = 0.0;

  	for (i = 0 ; i < 3 ; i++)
      	dist += pow((coord1[i] - coord2[i]), 2);

	return sqrt(dist);
}
