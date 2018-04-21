#include "functions.h"
#include "cRMSD.h"

// Function that returns True if a specific conformation is centroid
int is_centroid(int *centroids, int index, int size)
{
    int i;
    
    for (i = 0 ; i < size ; i++)
        if (centroids[i] == index)
            return 1;
        
    return 0;
}


// Function that generates a random number following Gaussian Distribution
double rand_gaussian()
{
	double r, y1, y2;

	do
	{
		y1 = (rand() / (RAND_MAX + 1.0)) * 2.0;
		y2 = (rand() / (RAND_MAX + 1.0)) * 2.0;  
		r = y1*y1 + y2*y2;	
	}
	while (r >= 1);

	return y1 * pow(((pow(r,-2) - 1)/r) , 1/2);
}


// Function that calculates the number of centroids that haven't been transposed ( or remotely transposed ) after an update 
int centroids_transposition(conformation **confs_array, int *centroids, int *centroids_old, int k_cluster, double **distances_array)
{
	int i, counter = 0;
	
	// for every single cluster
	for(i = 0 ; i < k_cluster ; i++)
		if (distances_array[centroids[i]][centroids_old[i]] < THRESHOLD)
			counter++;
			
	return counter;
}


// Auxiliary function that creates an array to store distances among every conformation
double **create_distances_array(conformation **confs_array, int num_conf, int N, int frechet)
{
	int i, j, y;
	double distance;
	double **array;
	
	array = (double**) malloc(num_conf * sizeof(double*));
	
	if (array == NULL)
	{	
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < num_conf ; i++)
	{
		array[i] = (double*) malloc(num_conf * sizeof(double));
			
		if (array[i] == NULL)
		{	
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}
	
	// Calculate distances among conformations
	for (i = 0 ; i < num_conf ; i++)
	{
		// Set main diagonial of distance array to 0
		array[i][i] = 0.0;
		
		for (j = 0 ; j < i ; j++)
		{
			distance = cRMSD(confs_array[i]->coordinates, confs_array[j]->coordinates, N, frechet);
			
			// Fill in the array
			array[i][j] = distance;
			array[j][i] = distance;
		}
	}

	return array;
}