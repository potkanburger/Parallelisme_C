
/*  tokenring example using PVM 3.4 
    - uses sibling() to determine the nb of spawned tasks (xpvm and pvm> ok)
    - uses group for token ring communication
*/

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include "pvm3.h"

#define GRPNAME "tokenring"

void matrix_load_paral(int me, int NPROC, int tids[], char nom[], double *tab, int N, int k);
void matrix_save_paral(int me, int NPROC, int tids[], char nom[], double *tab, int N, int k);
void matrix_display ( double *tab,int  N, int M );
void gauss ( double *tab, int N );

main(int argc, char ** argv) {

/*--------------------------------------------------------------------------*/
/*           										variables								                    */

    int NPROC = 2;							/* default nb of proc */
    int mytid;                  /* my task id */
    int *tids;                  /* array of task id */
    int me;                     /* my process number */
    int i;

/*--------------------------------------------------------------------------*/
/*           								PVM configuration					                      */

    /* enroll in pvm */
    mytid = pvm_mytid();

    /* determine the size of my sibling list */
    NPROC = pvm_siblings(&tids); 
    /* WARNING: tids are in order of spawning, which is different from
       the task index JOINING the group */

    me = pvm_joingroup( GRPNAME ); /* me: task index in the group */
    pvm_barrier( GRPNAME, NPROC );
    pvm_freezegroup ( GRPNAME, NPROC );
    for ( i = 0; i < NPROC; i++) tids[i] = pvm_gettid ( GRPNAME, i); 

/*--------------------------------------------------------------------------*/
/*           								starting the algo      					                */
     dowork( me, NPROC, tids, argc, argv );

/*--------------------------------------------------------------------------*/
/*           								exiting PVM group                  					    */
     pvm_lvgroup( GRPNAME );
     pvm_exit();
}


dowork( int me, int NPROC, int tids[], int argc, char ** argv) {

/*--------------------------------------------------------------------------*/
/*           						 variables Algorithm						                    */

		int N; // taille matrice 
		int i, j, k;
		double *tab, pivot;
		char nom[255];
		FILE *f;
		struct timeval tv1, tv2;	/* for timing */
		int duree;


		//Check args
		if (argc != 3){
			printf ("Usage: %s <matrix size> <matrix name>\n", argv[0]);
	
			exit (-1);
		} 

		//Malloc table for matrice
		N = atoi ( argv[1] );
		if ( (tab=malloc(N*(N/NPROC)*sizeof(double))) == NULL ) {
			printf ("Cant malloc %d bytes\n", N*(N/NPROC)*sizeof(double));
			exit (-1);
		}


		//recup nom de la matrice
		strcpy(nom, argv[2]);

		
		//gettimeofday( &tv1, (struct timezone*)0 );
		matrix_load_paral (me, NPROC, tids, nom , tab, N, 2 );
		//gettimeofday( &tv2, (struct timezone*)0 );
		//duree = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
		//printf ("loading time: %10.8f sec.\n", duree/1000000.0 );

		//gettimeofday( &tv1, (struct timezone*)0 );
		//gauss ( tab, N );
		//gettimeofday( &tv2, (struct timezone*)0 );
		//duree = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
		//printf ("computation time: %10.8f sec.\n", duree/1000000.0 );

		//sprintf ( nom+strlen(nom), ".result" );
		

		matrix_save_paral (me, NPROC, tids, nom, tab, N, 2 );

		//printf ("proc %i : \n", me);
		//matrix_display(tab, N/NPROC, N);
		//fflush(stdout);

}

void matrix_load_paral (int me, int NPROC, int tids[], char nom[], double *tab, int N, int k ) {

	double *token;
	double *line, *first_line;
	int src, dest;
	int count  = N-k;
	int stride = 1;
	double *temp;

  FILE *f;
  int i,j,p;
	int l0 = 0;

	//[PROC 0]
	if(me==0){
		double *buffer;

		if ( (token=malloc((N-k)*sizeof(double))) == NULL ) {
			printf ("Cant malloc %d bytes\n", (N-k)*(N-k)*sizeof(double));
			exit (-1);
		}

		if ( (temp=malloc((N-k)*sizeof(double))) == NULL ) {
			printf ("Cant malloc %d bytes\n", (N-k)*(N-k)*sizeof(double));
			exit (-1);
		}

		if ((f = fopen (nom, "r")) == NULL) { perror ("matrix_load : fopen "); }

		for (i=0; i<N; i++) {
			int me_cible = (i-k) % NPROC;
			if(i >= k){
			

				//lire ligne
				for (j=0; j<N; j++) {
					if(j >= k){
						fscanf (f, "%lf", (token+j-k) );
					}
					else
					{
						fscanf (f, "%lf", temp );
					}
				}
			
				//si me_cible = 0 on stocke la ligne
				if(me_cible == 0){
					memcpy(tab+(N-k)*l0, token, sizeof(double)*(N-k));		
					l0 ++;
				}

				if(i == k)
				{
					for (p=0; p<NPROC; p++)
					{
						dest = tids[p];
						
					  pvm_initsend( PvmDataDefault );
					  pvm_pkdouble( tab, count, stride );
					  pvm_send( dest, 0 );
					}
				}

				//sinon on l'envoie
				else{
					dest = tids[me_cible];
			    pvm_initsend( PvmDataDefault );
			    pvm_pkdouble( token, count, stride );
			    pvm_send( dest, 0 );

				}

			}

			else
			{
				for (j=0; j<N; j++) {
					fscanf (f, "%lf", temp );
				}
			}	
		
		}

		printf ("Close %s \n", nom);
		fflush(stdout);
		fclose (f);

	}
	

	//[PROC !0]
	else{

		if ( (first_line=malloc((N-k)*sizeof(double))) == NULL ) {
			printf ("Cant malloc %d bytes\n", (N-k)*(N-k)*sizeof(double));
			exit (-1);
		}
		src = tids[0];
	  pvm_recv( src, -1);
		pvm_upkdouble( first_line, count, stride );
		printf ("P%i reçu first line :\n", me);
		matrix_display(first_line, 1, N-k);
		fflush(stdout);

		for (i=0; i<(N-k)/NPROC; i++) {
		  pvm_recv( src, -1);
			pvm_upkdouble( tab+(N-k)*i, count, stride );
			printf ("P%i reçu :\n", me);
			matrix_display(tab+(N-k)*i, 1, N-k);
			fflush(stdout);
		}
		i++;
		if(me <= ((N-k)%NPROC)-1)
		{
		  pvm_recv( src, -1);
			pvm_upkdouble( tab+(N-k)*i, count, stride );
			printf ("P%i reçu :\n", me);
			matrix_display(tab+(N-k)*i, 1, N-k);
			fflush(stdout);
		}

	}


}

// sauvegarder les ligne du fichier par entrelacement

void matrix_save_paral (int me, int NPROC, int tids[], char nom[], double *tab, int N, int k ) {

	double *token;
	int src, dest, me_called;
	int count  = N-k;
	int stride = 1;

  FILE *f;
  int i,j;
	int l0 = 0;
	char c;

	if ( (token=malloc((N-k)*sizeof(double))) == NULL ) {
		printf ("Cant malloc %d bytes\n", (N-k)*(N-k)*sizeof(double));
		exit (-1);
	}

	//[PROC 0]
	if(me==0){
		
		strcat(nom,".result");

		f = fopen(nom, "a+");
		rewind(f);

		if (!f)
    	perror ("cant create file");

		printf ("created %s\n", nom);
		fflush(stdout);

		for (i=0; i<N; i++) {



			if(i>=k)
			{
				me_called = (i-k)%NPROC;

				printf ("P0 : i = %i, me_called = %i\n", i, me_called);
				fflush(stdout);

				if(me_called != 0 )
				{
					src = tids[me_called];

					printf ("P0 : en attente de me_%i\n", me_called);
					fflush(stdout);

					pvm_recv( src, 0 );
					pvm_upkdouble( token, count, stride );

					printf ("P0 : reçu de me_%i\n", me_called);
					fflush(stdout);


					printf ("ecrit ligne %i : \n", i);
					matrix_display(token, 1, N-k);
					fflush(stdout);
					for (j=0; j<N-k; j++) {
						if(k >=k)
							{fprintf (f, "%8.2f ", *(token+j));}
						else
							{fscanf (f, "%lf", NULL );}
						
					}
				}
				else
				{
					printf ("ecrit ligne %i : \n", i);
					matrix_display((tab+l0*N), 1, N-k);
					fflush(stdout);
					for (j=0; j<N-k; j++) {
						if(k >=k)
							{fprintf (f, "%8.2f ", *(tab+l0*N+j));}
						else
							{fscanf (f, "%lf", NULL );}
					}
					l0++;
				}

		  	fprintf (f, "\n");


			}
			else
			{
				printf ("saute ligne %i\n", i);
				fflush(stdout);
				while(c != '\n')
				{
					fscanf (f, "c", c );
					printf ("lut : %c \n", c);
					fflush(stdout);
				}
			}

		}
		fclose (f);

		printf ("fin \n");
		fflush(stdout);
	}
  
	//[PROC ! 0]
	else{
		for (i=0; i<(N-k)/NPROC; i++) {
			dest = tids[0];
		  token = tab+(N-k)*i;
			
			printf ("sending from %i :", me);
			matrix_display(token, 1, N-k);
			fflush(stdout);

		  pvm_initsend( PvmDataDefault );
		  pvm_pkdouble( token, count, stride );
		  pvm_send( dest, 0 );

			printf ("sent\n");
			fflush(stdout);
		}
		i++;
		if(me <= ((N-k)%NPROC)-1)
		{
			dest = tids[0];
		  token = tab+(N-k)*i;
			
			printf ("sent from %i :", me);
			matrix_display(token, 1, N-k);
			fflush(stdout);

		  pvm_initsend( PvmDataDefault );
		  pvm_pkdouble( token, count, stride );
		  pvm_send( dest, 0 );
		}
	}

}

void gauss ( double* tab, int N ){
		
}

void matrix_display ( double *tab,int  N, int M ) {
  int i,j;

  for (i=0; i<N; i++) {
    for (j=0; j<M; j++) {
      printf ("%8.2f ", *(tab+i*N+j) );
    }
    printf ("\n");
  }

}


