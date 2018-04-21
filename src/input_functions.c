#include "functions.h"
#include "input_functions.h"

// Function that reads the data from the input file and stores them into an array
conformation **parse_input(char *input_file)
{
	int i, j, numConform, N;
	char line[LINE_LENGTH];
	char *data;
	
	conformation **confs_array;
	
	// Open input file
	FILE *input = fopen(input_file, "r");
	
	// Check if input file has opened successfully
  	if (input == NULL)
  	{
    	printf("Fopen: error opening input file!\n");
      	exit(1);
    }
	
	// Get numConform from input file
  	fgets(line, sizeof(line), input);
	data = strtok(line,"\n"); 
  	numConform = atoi(data);
	
	// Get N from input file
	fgets(line, sizeof(line), input);
	data = strtok(line,"\n");
	N = atoi(data);
	
	// Array that stores dataset's conformations
	confs_array = (conformation**) malloc(numConform * sizeof(conformation*));
	
	if (confs_array == NULL)
	{
		printf("Malloc: memory allocation error (input_functions.c)!\n");
		exit(3);
	}
	
	// For every single conformation
	for(i = 0; i < numConform; i++)
	{
		
		conformation *new_conf = (conformation*) malloc(sizeof(conformation));
	
		if (new_conf == NULL)
		{
			printf("Malloc: memory allocation error (input_functions.c)!\n");
			exit(3);
		}

		// Assign conformation's number of amino acids
		new_conf->N = N;
		// Assign conformation's ID
		new_conf->ID = i;
		// Allocate memory for array to store the points
		new_conf->coordinates = (double**) malloc((new_conf->N) * sizeof(double*));
		
		if (new_conf->coordinates == NULL)
		{
			printf("Malloc: memory allocation error (input_functions.c)!\n");
			exit(3);
		}

		// Parse every single amino acid of conformation
		for (j = 0 ; j < N ; j++)
		{	
			// Allocate memory for array to store each coordinate for every point (related to dimension {2,3,4})
			new_conf->coordinates[j] = (double*) malloc(3 * sizeof(double));
					
			if (new_conf->coordinates[j] == NULL)
			{
				printf("Malloc: memory allocation error (input_functions.c)!\n");
				exit(3);
			}

			fgets(line, sizeof(line), input);
			
			// Get x coordinate
			data = strtok(line,"\t");
			new_conf->coordinates[j][0] = atof(data);
			
			// Get y coordinate
			data = strtok(NULL,"\t");
			new_conf->coordinates[j][1] = atof(data);
			
			// Get z coordinate
			data = strtok(NULL,"\n");
			new_conf->coordinates[j][2] = atof(data);
		}
		
		// Store the new conformation struct in array
		confs_array[i] = new_conf;
	}

	// Close input file
	fclose(input);
	
	return confs_array;
}


// Auxiliary function to get the number of conformations
int get_conf(char *input_file)
{

	char line[LINE_LENGTH];
	char *data;

	// Open input file
	FILE *input = fopen(input_file, "r");
	
	// Check if input file has opened successfully
  	if (input == NULL)
  	{
    	printf("Fopen: error opening input file!\n");
      	exit(1);
    }
	
	// Get numConform from input file
  	fgets(line, sizeof(line), input);
	data = strtok(line,"\n"); 

	// Close input file
	fclose(input);
	
	return atoi(data);
}

// Auxiliary function to get the number of points, that a conformation consists of
int get_N(char *input_file)
{

	char line[LINE_LENGTH];
	char *data;

	// Open input file
	FILE *input = fopen(input_file, "r");
	
	// Check if input file has opened successfully
  	if (input == NULL)
  	{
    	printf("Fopen: error opening input file!\n");
      	exit(1);
    }
	
	// Skip first line
  	fgets(line, sizeof(line), input);
	
	// Get N from input file
	fgets(line, sizeof(line), input);
	data = strtok(line,"\n"); 

	// Close input file
	fclose(input);
	
	return atoi(data);
}