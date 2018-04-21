#include "functions.h"
#include "silhouette.h"

// Function that calculates the silhouette for each cluster and the total silhouette of the whole clustering procedure
double silhouette(conformation **confs_array, int *centroids, cluster **clusters, int k_cluster, double **distances_array)
{
	int i, j, k, index;
	double s, sum, avg_a, avg_b, s_total, min_distance, sil;
	double *total_s ;
	
	
	total_s = (double*) malloc((k_cluster+1) * sizeof(double));
	
	if (total_s == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// For every cluster..
	for (i = 0 ; i < k_cluster ; i++)
	{
		s_total = 0;
		
		// For every object in this cluster..
		for (j = 0 ; j < clusters[i]->cluster_size ; j++)
		{
			// Calculate the distance between this object and its centroid
			sum = distances_array[centroids[i]][clusters[i]->confs_in_cluster[j]];
            			
			// Calculate the sum distance among this object and every other object in this cluster			
			for (k = 0 ; k < clusters[i]->cluster_size ; k++)
				sum += distances_array[clusters[i]->confs_in_cluster[j]][clusters[i]->confs_in_cluster[k]];
				
			
			// Calculate the average distance for this specific object
			avg_a = (double) sum / clusters[i]->cluster_size;

			
			min_distance = INFINITY;
			
			// Calculate the second nearest centroid from this specific object
			for (k = 0 ; k < k_cluster ; k++)
				if (k != i)
				{	
					sum = distances_array[clusters[i]->confs_in_cluster[j]][centroids[k]];
					
					if (sum < min_distance)
					{
						min_distance = sum;
						index = k;
					}
				}

			sum = min_distance;
			
			// And for this specific object calculate its shilhoutte for the second nearest centroid
			for (k = 0 ; k < clusters[index]->cluster_size ; k++)
				sum += distances_array[clusters[i]->confs_in_cluster[j]][clusters[index]->confs_in_cluster[k]];
				
						
			// Calculate the average distance for this specific object
			avg_b = (double) sum / (clusters[index]->cluster_size);
			
			if (avg_a >= avg_b)
				s = (double) (avg_b - avg_a) / avg_a;
			else
				s = (double) (avg_b - avg_a) / avg_b;
			
			s_total += s;		
		}
	
		// Update total_s array with #k_clusters shilhouttes ( for every cluster respectively )
		total_s[i] = s_total/clusters[i]->cluster_size;	
	}
	
	total_s[k_cluster] = 0;
	for (i = 0 ; i < k_cluster ; i++)
		total_s[k_cluster] += total_s[i];
	

	// Update total_s array's last position that is the mean shilhoutte of all clustering procedure
	sil = (double) (total_s[k_cluster] / k_cluster);

	free(total_s);
	
	return sil;
}