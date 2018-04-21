typedef struct _conformation_
{
	int N;
	int ID;
	double **coordinates;
} conformation;


typedef struct _cluster_
{
    int cluster_size;
    int *confs_in_cluster;
}cluster;