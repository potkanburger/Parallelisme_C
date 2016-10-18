#include "mpi.h"
#include <stdio.h>
int main( int argc, char *argv[] )
{
	int rank, size;
	char message[20];
	MPI_Status status;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	printf( "I am %d of %d\n", rank, size );
	//
	if (rank == 0) {
		/* code for process zero */
		strcpy(message,"Hello, there");
		
		/*MPI_Send(message, strlen(message), MPI_CHAR, 1, 99,
		MPI_COMM_WORLD);*/
		
	} 
	MPI_Bcast(&message, sizeof(char)*sizeof(message), MPI_CHAR, 0, MPI_COMM_WORLD);
	
	if(rank>0){
	/* code for other processes */
		//MPI_Recv(message, 20, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		printf("%d received :%s:\n", rank, message);
	}

	MPI_Finalize();
	return 0;
}


