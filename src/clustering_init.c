#include "functions.h"
#include "clustering_init.h"

// Implementation of K-means++
int *Init_1(conformation **confs_array, int k_cluster, int num_conf, double **distances_array) {
	int i, j, z, index, max_P;
    int *centroids;
	double distance, x, max;
	double *D, *P;

	// Array with min distances
	D = (double*) malloc(num_conf * sizeof(double));

	if (D == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	// Array with probabilities
	P = (double*) malloc(num_conf * sizeof(double));

	if (P == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	centroids = (int*) malloc(k_cluster * sizeof(int));

	if (centroids == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	//Choose first centroid uniformly at random
	index = (int) (rand_gaussian() * num_conf);

	//Assign it into array with centroids
	centroids[0] = index;

	// Initialization of arrays
	for(i = 1; i < k_cluster ; i++) {
		for(j = 0; j < num_conf; j++) {
			D[j] = INFINITY;
			P[j] = 0.0;
		}

		// For every single curve
		for(j = 0; j < num_conf ; j++) {
			// If its not already a centroid
			if (!is_centroid(centroids, j, i)) {
				// Calculate the distances from every centroid
				for(z = 0 ; z < i ; z++) {
					distance = distances_array[centroids[z]][j];

                    // and keep the minimum distance
                    if (distance < D[j])
                        D[j] = distance;
				}
			}
		}

		max = -INFINITY;

		// Find the max value of minimum distances in D
		for(j = 0; j < num_conf; j++)
			if (D[j] != INFINITY && D[j] > max)
				max = D[j];

		// Creation of P(r) array
		for(j = 0; j < num_conf; j++) {
			if (D[j] != INFINITY) {
				for(z = 0; z < j; z++) {
					if (D[z] != INFINITY)
						P[j] += pow(D[z]/max, 2);
				}
			}
		}

		// Keep in max_P variable the max value of P array
		for(j = num_conf-1 ; j >= 0; j--)
			if (P[j] != 0) {
				max_P = P[j];
				break;
			}
		// Choose randomly a real number in space [0, max_P]
		x = ((float)rand()/(float)(RAND_MAX)) * max_P;

		for(j = 0; j < num_conf; j++) {				
			if (x <= P[j]){
				//this is the curve that we'll choose for centroid
				centroids[i] = j;
				break;
			}
		}
	}

	free(D);
	free(P);

	return centroids;
}
