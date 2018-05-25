#include "functions.h"
#include "clustering_update.h"

// Implementation of Partitioning Around Medoids (PAM) – Improved update
void PAM(conformation **confs_array, int *centroids, cluster **clusters, int k_cluster, double **distances_array) {
	int i, j, z, min_index;
	double min;
	double *distances;


	// For every cluster..
	for (i = 0 ; i < k_cluster ; i++) {
		distances = (double*) malloc(clusters[i]->cluster_size * sizeof(double));

		if (distances == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		// For every object in this cluster..
		for (j = 0 ; j < clusters[i]->cluster_size ; j++) {
			distances[j] = 0;

			// Find the sum distance with every other object, and update the array distances..
			for (z = 0 ; z < clusters[i]->cluster_size ; z++)
				distances[j] += distances_array[clusters[i]->confs_in_cluster[j]][clusters[i]->confs_in_cluster[z]];

		}

		min = INFINITY;

		// Find out the min sum distance and define this object as the new centroid
		for (j = 0 ; j < clusters[i]->cluster_size ; j++) {
			if (distances[j] < min) {
				min = distances[j];
				min_index = j;
			}
		}

		// Update the array with centroids
		centroids[i] = clusters[i]->confs_in_cluster[min_index];

		free(distances);
	}
	return;
}
