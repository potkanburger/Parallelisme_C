#include "mpi.h"
#include <stdio.h>
#include <math.h>

int main( int argc, char *argv[] )
{
	int rank, size;
	MPI_Status status;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );

	int maxn = 32;
	double err = 0.0;
	double convergence = 1.0;
	double convergence_locale = 0.0;

	double dsq_size = sqrt(size);
	if ((int)dsq_size != dsq_size){
		return 0;
	}

	int sq_size = (int)dsq_size;

	if (maxn%sq_size != 0){
		maxn += sq_size - maxn%sq_size;
	}
	

	double f[(maxn/sq_size)+2][(maxn/sq_size)+2];
	double fnew[(maxn/sq_size)+2][(maxn/sq_size)+2];
	
	double xrecuG[(maxn/sq_size)+2];
	double xrecuD[(maxn/sq_size)+2];
	double xrecuH[(maxn/sq_size)+2];
	double xrecuB[(maxn/sq_size)+2];

	double tmpColonne[(maxn/sq_size)+2];

	int i,j;
	
	for (i = 0; i<(maxn/sq_size)+2; i++)
	{
		f[0][i] = -1;
		f[(maxn/sq_size)+1][i] = -1;
		f[i][0] = -1;
		f[i][(maxn/sq_size)+1] = -1;
	}

	for (i=1; i<(maxn /sq_size)+1;i++)
	{
		for (j=1;j<(maxn/sq_size)+1;j++)
		{
			f[i][j] = rank;
		}
	}
	
	
	// Si pas sur la première ligne, on envoie à la ligne du dessus et recoit d'elle
	if (rank >= sq_size){
		MPI_Send(&f[1], (maxn/sq_size)+2, MPI_DOUBLE, rank-sq_size, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuH, (maxn/sq_size)+2, MPI_DOUBLE, rank-sq_size, 99, MPI_COMM_WORLD, &status);
	}

	// Si pas sur la dernière ligne, on envoie à la ligne du dessous et recoit d'elle
	if (rank < size-sq_size){
		MPI_Send(&f[(maxn/sq_size)], (maxn/sq_size)+2, MPI_DOUBLE, rank+sq_size, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuB, (maxn/sq_size)+2, MPI_DOUBLE, rank+sq_size, 99, MPI_COMM_WORLD, &status);
	}

	// Si pas sur la première colonne, on envoie à la colonne de gauche et recoit d'elle
	if (rank%sq_size != 0){

		for (i=0; i<(maxn/sq_size)+2; i++){
			tmpColonne[i] = f[i][1];
		}

		MPI_Send(&tmpColonne, (maxn/sq_size)+2, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuG, (maxn/sq_size)+2, MPI_DOUBLE, rank-1, 99, MPI_COMM_WORLD, &status);
	}

	// Si pas sur la dernière colonne, on envoie à la colonne de droite et recoit d'elle
	if (rank%sq_size != sq_size-1){

		for (i = 0; i<(maxn / sq_size) + 2; i++){
			tmpColonne[i] = f[i][(maxn/sq_size)];
		}

		MPI_Send(&tmpColonne, (maxn/sq_size)+2, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD);
		MPI_Recv(&xrecuD, (maxn/sq_size)+2, MPI_DOUBLE, rank+1, 99, MPI_COMM_WORLD, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);


	// Si pas sur la première ligne, on recupere la ligne du dessus
	if (rank >= sq_size){
		for (i=0; i<(maxn/sq_size)+2;i++)
		{
			f[0][i] = xrecuH[i];
		}
	}

	// Si pas sur la dernière ligne, on recupere la ligne du dessous
	if (rank < size - sq_size){
		for (i=0; i<(maxn/sq_size)+2; i++)
		{
			f[(maxn/sq_size)+1][i] = xrecuB[i];
		}
	}

	// Si pas sur la première colonne, on recupere la colonne de gauche
	if (rank%sq_size != 0){

		for (i=0;i<(maxn/sq_size)+2; i++){
			f[i][0] = xrecuG[i];
		}

	}

	// Si pas sur la dernière colonne, on recupere la colonne de droite
	if (rank%sq_size != sq_size - 1){

		for (i=0; i<(maxn/sq_size)+2;i++){
			f[i][(maxn/sq_size)+1] = xrecuD[i];
		}
	}


		
	for(i=0;i<(maxn/sq_size)+2;i++){
		for (j=0; j<(maxn/sq_size)+2; j++){
			fnew[i][j] = f[i][j];
		}
	}
	
	
	
	
	while(convergence>0.001) {
				
		for(i=1;i<(maxn/sq_size)+1;i++){
			for (j=1;j<(maxn/sq_size)+1;j++){
				fnew[i][j] = 0.25*(f[i+1][j]+f[i-1][j]+f[i][j+1]+f[i][j-1]);
			}
		}
		
		err = 0.0;
		for(i=1;i<(maxn/sq_size)+1;i++){
			for (j=1;j<(maxn/sq_size)+1;j++){
				err += (fnew[i][j]-f[i][j])*(fnew[i][j]-f[i][j]);	
			}
		}

		convergence_locale = sqrt(err);
		
		MPI_Allreduce(&convergence_locale, &convergence, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
			
		for(i=0;i<(maxn/sq_size)+2;i++){
			for(j=0; j<(maxn/sq_size)+2; j++){
				f[i][j] = fnew[i][j];
			}
		}



		// Si pas sur la première ligne, on envoie à la ligne du dessus et recoit d'elle
		if(rank >= sq_size){
			MPI_Send(&fnew[1], (maxn / sq_size) + 2, MPI_DOUBLE, rank - sq_size, 99, MPI_COMM_WORLD);
			MPI_Recv(&xrecuH, (maxn / sq_size) + 2, MPI_DOUBLE, rank - sq_size, 99, MPI_COMM_WORLD, &status);
		}

		// Si pas sur la dernière ligne, on envoie à la ligne du dessous et recoit d'elle
		if(rank < size - sq_size){
			MPI_Send(&fnew[(maxn / sq_size)], (maxn / sq_size) + 2, MPI_DOUBLE, rank + sq_size, 99, MPI_COMM_WORLD);
			MPI_Recv(&xrecuB, (maxn / sq_size) + 2, MPI_DOUBLE, rank + sq_size, 99, MPI_COMM_WORLD, &status);
		}

		// Si pas sur la première colonne, on envoie à la colonne de gauche et recoit d'elle
		if(rank%sq_size != 0){

			for(i = 0; i < (maxn / sq_size) + 2; i++){
				tmpColonne[i] = fnew[i][1];
			}

			MPI_Send(&tmpColonne, (maxn / sq_size) + 2, MPI_DOUBLE, rank - 1, 99, MPI_COMM_WORLD);
			MPI_Recv(&xrecuG, (maxn / sq_size) + 2, MPI_DOUBLE, rank - 1, 99, MPI_COMM_WORLD, &status);
		}

		// Si pas sur la dernière colonne, on envoie à la colonne de droite et recoit d'elle
		if(rank%sq_size != sq_size - 1){

			for(i = 0; i < (maxn / sq_size) + 2; i++){
				tmpColonne[i] = fnew[i][(maxn / sq_size)];
			}

			MPI_Send(&tmpColonne, (maxn / sq_size) + 2, MPI_DOUBLE, rank + 1, 99, MPI_COMM_WORLD);
			MPI_Recv(&xrecuD, (maxn / sq_size) + 2, MPI_DOUBLE, rank + 1, 99, MPI_COMM_WORLD, &status);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// Si pas sur la première ligne, on recupere la ligne du dessus
		if(rank >= sq_size){
			for(i = 0; i < (maxn / sq_size) + 2; i++)
			{
				f[0][i] = xrecuH[i];
			}
		}

		// Si pas sur la dernière ligne, on recupere la ligne du dessous
		if(rank < size - sq_size){
			for(i = 0; i < (maxn / sq_size) + 2; i++)
			{
				f[(maxn / sq_size) + 1][i] = xrecuB[i];
			}
		}

		// Si pas sur la première colonne, on recupere la colonne de gauche
		if(rank%sq_size != 0){

			for(i = 0; i < (maxn / sq_size) + 2; i++){
				f[i][0] = xrecuG[i];
			}

		}

		// Si pas sur la dernière colonne, on recupere la colonne de droite
		if(rank%sq_size != sq_size - 1){

			for(i = 0; i < (maxn / sq_size) + 2; i++){
				f[i][(maxn / sq_size) + 1] = xrecuD[i];
			}
		}

	}


	if (rank > 0){
		MPI_Send(&f, ((maxn / sq_size) + 2)*((maxn / sq_size) + 2), MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
	}


	if (rank == 0){
		FILE * file;
		int k;
		int result[maxn][maxn];
		int ligne;
		int colonne;

		file = fopen("result_decomposition2D", "w+");
		double recu[(maxn / sq_size) + 2][(maxn / sq_size) + 2];
		double bord = -1.00;


		for (i=1; i<(maxn/sq_size)+1); i++){
			for (j = 1; j < (maxn / sq_size) + 1; j++){
				result[i - 1][j - 1] = f[i][j];
			}
		}

		for(k=1;k<size;k++){
			MPI_Recv(&recu, ((maxn/sq_size)+2)*((maxn/sq_size)+2), MPI_DOUBLE, k, 99, MPI_COMM_WORLD, &status);
			
			ligne = (k/sq_size)*(maxn/sq_size);
			colonne = (k%sq_size)*(maxn/sq_size);

			for (i = 0; i<maxn/sq_size; i++){
				for (j = 0; j < maxn/sq_size; j++){
					result[ligne+i][colonne+j] = recu[i+1][j+1];
				}
			}
		}

		for (i=0; i<maxn; i++)
		{
			for (j = 0; j<maxn; j++)
			{
				fprintf(file, "%.2f ", result[i][j]);
			}
			fprintf(file, "\n");
		}
		fflush(stdout);
		fclose(file);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	return 0;
}
