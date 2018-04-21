#include "functions.h"
#include "cRMSD.h"
#include "metric_functions.h"

// Implementation of c-RMSD algorithm
double cRMSD(double **coord1, double **coord2, int N, int frechet)
{
	int i, j;
	int *ipiv;
	double det, F_norm, distance;
	double *c, *singular, *V_transpose, *U, *Q, *superb, *tmp_matrix, *tmp1_matrix, *tmp2_matrix;
	double **centroids, **coord1_translate, **coord2_translate, **coord1_transpose;

	coord1_translate = (double**) malloc(N * sizeof(double*));
	
	if (coord1_translate == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	coord2_translate = (double**) malloc(N * sizeof(double*));
	
	if (coord2_translate == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}

	for (i = 0 ; i < N ; i++)
	{
		coord1_translate[i] = (double*) malloc (3 * sizeof(double));
		
		if (coord1_translate[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	
		coord2_translate[i] = (double*) malloc (3 * sizeof(double));
		
		if (coord2_translate[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}

	// Find the centroid of each pointset
	centroids = find_centroids(coord1, coord2, N);
	
	for (i = 0 ; i < N ; i++)
	{
		for (j = 0 ; j < 3 ; j++)
		{
			// X translated
			coord1_translate[i][j] = coord1[i][j] - centroids[0][j];
			// Y translated
			coord2_translate[i][j] = coord2[i][j] - centroids[1][j];
		}
	}
	
	// Create tranposed matrix of first pointset
	coord1_transpose = transpose(coord1_translate, N);
	
	
	// Multiply transposed X with Y
	// C matrix : result of X,Y matrices multiplication
	c = (double*) malloc(3 * 3 * sizeof(double));
	
	if (c == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// Mutiply X transpose and Y matrices
	tmp1_matrix = create_1D_array(coord1_transpose, 3, N);
	tmp2_matrix = create_1D_array(coord2_translate, N, 3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans , CblasNoTrans, 3, 3, N, 1.0, tmp1_matrix, 
				N, tmp2_matrix, 3, 0.0, c, 3);

	free(tmp1_matrix);
	free(tmp2_matrix);
	
	// SVD
	singular = LAPACKE_malloc(3 * sizeof(double));
	
	if (singular == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	V_transpose = LAPACKE_malloc(3 * 3 * sizeof(double));	
	
	if (V_transpose == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	U = LAPACKE_malloc(3 * 3 * sizeof(double));
	
	if (U == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	superb = LAPACKE_malloc(2 * sizeof(double));
	
	if (superb == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// SVD function from Lapacke
	LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', 
					3, 3, c, 3, singular, U, 
					3, V_transpose, 3, superb);

	// If σ3 is negative , we have to stop the procedure
	if (singular[2] < 0 )
	{
		printf("Degenerated Data...procedure failed!\n");
		return -1.0;
	}
	
	// Create Q matrix out of U and V_transpose multiplication
	Q = LAPACKE_malloc(3 * 3 * sizeof(double));	

	if (Q == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	// Multiplication of U and V transpose matrices
	cblas_dgemm(CblasRowMajor, CblasNoTrans , CblasNoTrans,
				3, 3, 3, 1.0, U, 
				3, V_transpose, 3, 0.0, Q, 3);


	// Calculate and check determinant of Q matrix
	ipiv = (int*) malloc(3 * sizeof(int));
		
	if (ipiv == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}	
		
	if(LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, Q, 3, ipiv) != 0)
	{
		printf("LAPACKE_dgetrf error!\n");
		return -2.0;
	}
		
	det = 1.0;
	for (i = 0; i < 3; i++)
	{
		 if (ipiv[i] != (i+1))
			 det = -det * Q[i*3+i];
		 else
			 det = det * Q[i*3+i];
	}
	
	//If det < 0
	if (det < 0) 
	{
		// The last col of U has to be multiplied with -1
		U[2] *= -1;
		U[5] *= -1;
		U[8] *= -1;
		
		// [U1 ,U2 ,-U3] * V transpose
		cblas_dgemm(CblasRowMajor, CblasNoTrans , CblasNoTrans,
				3, 3, 3, 1.0, U, 
				3, V_transpose, 3, 0.0, Q, 3);
	}

	
	tmp_matrix = LAPACKE_malloc(N * 3 * sizeof(double));
	
	if (tmp_matrix == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}	
		
	// tmp_matrix = X * Q
	tmp1_matrix = create_1D_array(coord1_translate, N, 3);
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans , CblasNoTrans,
				N, 3, 3, 1.0, tmp1_matrix, 
				3, Q, 3, 0.0, tmp_matrix, 3);
	
	free(tmp1_matrix);

	// Use either Frechet or Frobenius norm to calculate the distance
	if (frechet)
	{			
		for (i = 0 ; i < N ; i++)
			for (j = 0 ; j < 3 ; j++)
				coord1_translate[i][j] = tmp1_matrix[i*3+j];
			
		distance = frechet_distance(coord1_translate, coord2_translate, N);
	}
	else
	{
		// tmp_matrix = tmp_matrix - (Y TRANSLATE)
		for (i = 0 ; i < N * 3 ; i++)
			tmp_matrix[i] -= coord2_translate[i/3][i%3];
			
		// Calculate Frobenius norm
		F_norm = LAPACKE_dlange(CblasRowMajor, 'F', N, 3, tmp_matrix, 3);
		distance = F_norm / sqrt(N);
	}
	
	
	// Free allocated memory of the whole content
	// ==========================================
	
	// Free memory of coord1_transpose
	for (i = 0 ; i < 3 ; i++)
		free(coord1_transpose[i]);	
	free(coord1_transpose);
	
	// Free memory of coord1_translate, coord2_trans
	for (i = 0 ; i < N ; i++)
	{
		free(coord1_translate[i]);
		free(coord2_translate[i]);
	}
	free(coord1_translate);
	free(coord2_translate);
	
	// Free memory of centroids
	for (i = 0 ; i < 2 ; i++)
		free(centroids[i]);
	free(centroids);
	
	free(c);
	free(ipiv);
	
	LAPACKE_free(singular);
	LAPACKE_free(V_transpose);
	LAPACKE_free(U);
	LAPACKE_free(superb);
	LAPACKE_free(Q);
	LAPACKE_free(tmp_matrix);

	return distance;
}


// Auxiliary Function that calculates the transpose matrix
double **transpose(double **coord, int N)
{
	int i, j;
	double **trans;
	
	trans = (double**) malloc(3 * sizeof(double*));
	
	if (trans == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < 3 ; i++)
	{
		trans[i] = (double*) malloc(N * sizeof(double));
	
		if (trans[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
	}
	
	//Transpose Matrix calculation
	for (i = 0 ; i < 3 ; i++)
		for (j = 0 ; j < N ; j++)
			trans[i][j] = coord[j][i];
	
	return trans;
}


// Auxiliary Function that transforms a 2d-Array into 1d-array
double *create_1D_array(double **coord, int row, int col)
{
	int i, j;
	double *coord_1D;

	coord_1D = (double*) malloc(row * col * sizeof(double));
	
	if (coord_1D == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < row ; i++)
		for (j = 0 ; j < col ; j++)
			coord_1D[i*col+j] = coord[i][j];

	return coord_1D;
}


// Auxiliary Function that finds out the centroids of 2 pointsets 
double **find_centroids(double **coord1, double **coord2, int N)
{
	int i, j;
	double **centroids;

	centroids = (double**) malloc(2 * sizeof(double*));

	if (centroids == NULL)
	{
		printf("Malloc: memory allocation error!\n");
		exit(3);
	}
	
	for (i = 0 ; i < 2 ; i++)
	{
		centroids[i] = (double*) malloc (3 * sizeof(double));
		
		if (centroids[i] == NULL)
		{
			printf("Malloc: memory allocation error!\n");
			exit(3);
		}
		
		// Initialize centroids to x = 0, y = 0, z = 0
		for (j = 0 ; j < 3 ; j++)
			centroids[i][j] = 0.0;
	}

	// Find centroids by calculating the mean values of each pointset
	for (i = 0 ; i < N ; i++)
	{
		for (j = 0 ; j < 3 ; j++)
		{	
			centroids[0][j] += coord1[i][j];
			centroids[1][j] += coord2[i][j];	
		}
	}

	for (j = 0 ; j < 3 ; j++)
	{
		centroids[0][j] = centroids[0][j] / (double) N;
		centroids[1][j] = centroids[1][j] / (double) N;
	}

	return centroids;
}