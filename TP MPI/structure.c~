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
	double xlocal[(maxn/size)+2][maxn];
	double xrecuA[maxn];
	double xrecuB[maxn];
	int i,j;
	for(i=0;i<maxn;i++)
	{
		xlocal[0][i] = -1;
		xlocal[(maxn/size)+1][i] = -1;	
	}
	for(i=1;i<(maxn/size)+1;i++)
	{
		for(j=0;j<maxn;j++)
		{
			xlocal[i][j] = rank;
		}
	}
	
	
	if(rank==0){
		MPI_Send(&xlocal[(maxn/size)], maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuB, maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD, &status);	
	}
	else if(rank==(size-1)){
		MPI_Send(&xlocal[1], maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuA, maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD, &status);
	}
	else{
		MPI_Send(&xlocal[1], maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuA, maxn, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD, &status);
		
		MPI_Send(&xlocal[(maxn/size)], maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuB, maxn, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD, &status);	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	if(rank==0){
		for(i=0;i<maxn;i++)
		{
			xlocal[(maxn/size)+1][i] = xrecuB[i];	
		}
	}
	else if(rank==(size-1)){
		for(i=0;i<maxn;i++)
		{
			xlocal[0][i] = xrecuA[i];	
		}
	}
	else{
		for(i=0;i<maxn;i++)
		{
			xlocal[0][i] = xrecuA[i];	
			xlocal[(maxn/size)+1][i] = xrecuB[i];	
		}
	}
	
	
	if(rank==1){
		printf("rank: %d\n", rank);
		for(i=0;i<(maxn/size)+2;i++)
		{
			for(j=0;j<maxn;j++)
			{
				printf("%.1f ", xlocal[i][j]);
			}
			printf("\n");		
		}
		printf("\n\n");
		fflush(stdout);		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
}
