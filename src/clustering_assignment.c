#include "functions.h"
#include "clustering_assignment.h"

// Implementation of Lloyd’s assignment (simplest approach)
cluster **Assignment_1(conformation **confs_array, int k_cluster, int *centroids, int num_conf, double **distances_array) {
	int i, j, ID;
	double min_D, distance;

	cluster **clusters;

	// Initilization of k clusters
	clusters = (cluster**) malloc(k_cluster * sizeof(cluster*));

	if (clusters == NULL) {
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	for(i = 0; i < k_cluster ; i++) {
		clusters[i] = (cluster*) malloc(sizeof(cluster));

		if (clusters[i] == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		clusters[i]->cluster_size = 1;
		clusters[i]->confs_in_cluster = (int*) malloc(sizeof(int));

		if (clusters[i]->confs_in_cluster == NULL) {
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}

		clusters[i]->confs_in_cluster[0] = centroids[i];
	}

	// For every curve of dataset
	for(i = 0; i < num_conf ; i++) {
		// If its not a centroid
		if (!is_centroid(centroids, i, k_cluster)) {
			min_D = INFINITY;

			// For every centroid
			for(j = 0; j < k_cluster ; j++) {
				distance = distances_array[centroids[j]][i];

				// And keep the minimum distance
				if (distance < min_D) {
					min_D = distance;
					ID = j;
				}
			}

			clusters[ID]->cluster_size++;
			clusters[ID]->confs_in_cluster = (int*) realloc(clusters[ID]->confs_in_cluster, clusters[ID]->cluster_size * sizeof(int));

			if (clusters[ID]->confs_in_cluster == NULL) {
				printf("Malloc: memory allocation error!\n");
				exit(3);
			}

			clusters[ID]->confs_in_cluster[clusters[ID]->cluster_size-1] = i;
		}
	}
	return clusters;
}
