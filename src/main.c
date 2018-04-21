#include "functions.h"
#include "input_functions.h"
#include "clustering_init.h"
#include "clustering_assignment.h"
#include "clustering_update.h"
#include "output_functions.h"
#include "silhouette.h"

// execution command: ./proteins -i datasets/bio_small_input.dat -frechet [optional]

int main(int argc, char* argv[])
{
	int i, j, k_cluster, best_k_cluster = 0, num_conf, N, frechet = 0;
    int *centroids, *centroids_old;
	double sil, best_sil = -INFINITY;
	double **distances_array;
	char input_file[PATH_LENGTH] = "NULL";
	
	conformation **confs_array;
    cluster **clusters, **best_clusters;
		
		
    // Parse arguments from command line
    if (argc >= 1 && argc <= 4)									
    {
        for (i = 1 ; i <= (argc-1) ; i++)
        {
            if (!strcmp(argv[i], "-i")) // If flag is -i ..
            {
                stpcpy(input_file, argv[++i]); // Copy of argument's value to var
                continue;
            }
            if (!strcmp(argv[i], "-frechet")) // If flag is -frechet ..
            {
                frechet = 1; // Set this var True
                continue;
            }
            // If flag of any argument is wrong .. exit!
            printf("Wrong input arguments! \n");
            return -2;
        }
    }
    else	// If number of arguments is wrong .. exit!
    {
        printf("Wrong number of arguments! \n");
        return -3;	
    }

    // Parse input from user
    if (!strcmp(input_file, "NULL"))
    {
        printf("Please insert the path to input file: ");
        scanf("%s", input_file);
    } 

	// Seed random number generator 
	srand(time(NULL));
	
	// Parse data from input file
	confs_array = parse_input(input_file);
    
    num_conf = get_conf(input_file);
	
	N = get_N(input_file);
	
	printf("\tCreating distances's array\n\n");
	distances_array = create_distances_array(confs_array, num_conf, N, frechet);
	
	
	for (k_cluster = num_conf/4 ; k_cluster > 2 ; k_cluster/=2)
	{
		printf("---Perform clustering with %d clusters---\n", k_cluster);
		// Clustering
		centroids_old = (int*) malloc(k_cluster * sizeof(int));
		
		if (centroids_old == NULL)
		{
			printf("Malloc: memory allocation error/n");
			exit(3);
		}

		// Initialize centroids
		centroids = Init_1(confs_array, k_cluster, num_conf, distances_array);

		// Assign conformations to nearest cluster
		clusters = Assignment_1(confs_array, k_cluster, centroids, num_conf, distances_array);

		// Update centroids
		PAM(confs_array, centroids, clusters, k_cluster, distances_array);

		do
		{
			for (i = 0 ; i < k_cluster ; i++)			
				centroids_old[i] = centroids[i];

			// Assign conformations to nearest cluster
			clusters = Assignment_1(confs_array, k_cluster, centroids, num_conf, distances_array);

			// Update centroids
			PAM(confs_array, centroids, clusters, k_cluster,distances_array);
			
		}while(centroids_transposition(confs_array, centroids, centroids_old, k_cluster, distances_array) != k_cluster);
		
		// Calculate silhouette for current cluster
		sil = silhouette(confs_array, centroids, clusters, k_cluster, distances_array);
		
		printf("\tSilhouette: %f\n\n", sil);
		
		// Keep cluster with the best silhouette
		if (sil > best_sil)
		{
			// Free previous
			for(i = 0 ; i < best_k_cluster ; i++)
			{
				free(best_clusters[i]->confs_in_cluster);
				free(best_clusters[i]);
			}
			free(best_clusters);
			
			
			best_sil = sil;
			best_k_cluster = k_cluster;
			
			best_clusters = (cluster**) malloc(best_k_cluster * sizeof(cluster*));
			
			if (best_clusters == NULL)
			{
				printf("Malloc: memory allocation error/n");
				exit(3);
			}
			
			for (i = 0 ; i < best_k_cluster ; i++)
			{
				best_clusters[i] = (cluster*) malloc(sizeof(cluster));
				
				if (best_clusters[i] == NULL)
				{
					printf("Malloc: memory allocation error/n");
					exit(3);
				}
				
				best_clusters[i]->cluster_size = clusters[i]->cluster_size;
				
				best_clusters[i]->confs_in_cluster = (int*) malloc(best_clusters[i]->cluster_size *sizeof(int));
				
				if (best_clusters[i]->confs_in_cluster == NULL)
				{
					printf("Malloc: memory allocation error/n");
					exit(3);
				}
				
				for (j = 0 ; j < best_clusters[i]->cluster_size ; j++)
					best_clusters[i]->confs_in_cluster[j] = clusters[i]->confs_in_cluster[j];
			}
		}
		
		for(i = 0 ; i < k_cluster ; i++)
		{
			free(clusters[i]->confs_in_cluster);
			free(clusters[i]);
		}
		free(clusters);
		
		free(centroids);
		free(centroids_old);	
	}	
	
	printf("\tWriting results to file\n");
	
	// Write results to output file
	if (frechet)
		output("frechet.dat", best_clusters, best_k_cluster, best_sil);
	else
		output("crmsd.dat", best_clusters, best_k_cluster, best_sil);
	

	// Free allocated memory for distances array
	for (i = 0 ; i < num_conf ; i++)
		free(distances_array[i]);
	free(distances_array);
	
	// Free allocated memory for conformations array
	for (i = 0 ; i < num_conf ; i++)
	{
		for (j = 0 ; j < N ; j++)
			free(confs_array[i]->coordinates[j]);
		free(confs_array[i]->coordinates);
		free(confs_array[i]);
	}
	free(confs_array);
	
	
	return 0; 
}
