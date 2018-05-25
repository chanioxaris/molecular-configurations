#include "functions.h"
#include "output_functions.h"
#include "quicksort.h"

// Function that appends to output file
void output(char *output_file, cluster **clusters, int k_cluster, double silhouette) {
  	int i, j;

  	// Open the output_file to append new information
	FILE *output = fopen(output_file, "a");

  	// Check if the file opened successfully
  	if (output == NULL) {
		printf("Fopen: error opening output file!\n");
      	exit(1);
    }

  	// Let's start writing in file..
	fprintf(output, "k:\t%d\n", k_cluster);

	fprintf(output, "s:\t%f\n", silhouette);

  	for (i = 0 ; i < k_cluster ; i++) {
		// Sort	IDs before exporting to file
		quickSort(clusters[i]->confs_in_cluster, 0, clusters[i]->cluster_size -1);

		for (j = 0 ; j < clusters[i]->cluster_size ; j++)
			fprintf(output, "%d\t", clusters[i]->confs_in_cluster[j]);
		fprintf(output, "\n");
	}

  	fprintf(output, "\n\n");

  	// When the writing is finished, we close the file!
	fclose(output);

	return;
}
