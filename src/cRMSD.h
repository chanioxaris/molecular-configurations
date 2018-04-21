#ifndef _CRMSD_
#define _CRMSD_

// Declare functions of cRMSD.c file
double cRMSD(double **, double **, int, int);

double **transpose(double **, int);

double *create_1D_array(double **, int, int);

double **find_centroids(double **, double **, int);

#endif