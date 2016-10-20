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
	int decoupage = 10000000;
	int decoupageParProc = decoupage/size;
	//printf( "I am %d of %d\n", rank, size );
	//
	
	/*if (rank == 0) {
		//code for process zero
		strcpy(message,"Hello, there");
		double aire = 0;
		//MPI_Send(message, strlen(message), MPI_CHAR, 1, 99, MPI_COMM_WORLD);
		double partie;
		double x;
		for(partie = 0 ; partie < decoupage ; partie++)
		{
			x = partie/decoupage;
			aire = aire + ((1/(1+x*x)) * (1.0/decoupage));
			
		}
		aire = aire *4;
		printf("keshnée :%lf\n", aire);
		
		
		
	}*/
	double tmptime = MPI_Wtime();
	double comp_time; 
	
	int init = rank*decoupageParProc;
	double aire = 0;
	double partie;
	double x;
	double resultat;
	for(partie = init ; partie < decoupageParProc*(rank+1) ; partie++)
	{
		x = partie/decoupage;
		aire = aire + ((1/(1+x*x)) * (1.0/decoupage));		
	}
	//printf("rank %d : aire %lf\n", rank, aire);
	
	//MPI_Bcast(&message, sizeof(char)*sizeof(message), MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Reduce(&aire, &resultat, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	comp_time = MPI_Wtime()-tmptime;
	if(rank==0){
	/* code for other processes */
		//MPI_Recv(message, 20, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		resultat = resultat*4;
		printf("Comp time for proc. 0 = %5.2f seconds\n", comp_time);
		printf("resultat %lf\n", resultat);
	}

	MPI_Finalize();
	return 0;
}


