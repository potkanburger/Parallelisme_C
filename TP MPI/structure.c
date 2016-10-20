#include "mpi.h"
#include <stdio.h>
int main( int argc, char *argv[] )
{
	int rank, size;
	MPI_Status status;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	int maxn = 20;
	double x[maxn][maxn];
	double xlocal[maxn][(maxn/size)+2];
	double xrecuA[maxn];
	double xrecuB[maxn];
	int i,j;
	for(i=0;i<maxn;i++)
	{
		xlocal[i][0] = -1;
		xlocal[i][(maxn/size)+1] = -1;	
	}
	for(j=1;j<(maxn/size)+1;j++)
	{
		for(i=0;i<maxn;i++)
		{
			xlocal[i][j] = rank;
		}
	}
	
	if(rank==0){
		MPI_Send(&xlocal[(maxn/size)], maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuB, maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD, &status);
		for(i=0;i<maxn;i++)
		{
			xlocal[i][(maxn/size)+1] = xrecuB[i];	
		}
	}
	else if(rank==(size-1)){
		MPI_Send(&xlocal[0], maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuA, maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD, &status);
		for(i=0;i<maxn;i++)
		{
			xlocal[i][0] = xrecuA[i];	
		}
	}
	else{
		MPI_Send(&xlocal[0], maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD);
		MPI_Send(&xlocal[(maxn/size)], maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuA, maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&xrecuB, maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD, &status);
		for(i=0;i<maxn;i++)
		{
			xlocal[i][0] = xrecuA[i];	
			xlocal[i][(maxn/size)+1] = xrecuB[i];	
		}
	}
	printf("rank: %d", rank);
	for(j=0;j<(maxn/size)+2;j++)
	{
		for(i=0;i<maxn;i++)
		{
			printf("%lf", xlocal[i][j]);
		}
	}
	printf("\n");
	
}
