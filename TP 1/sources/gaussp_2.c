
/*  tokenring example using PVM 3.4 
    - uses sibling() to determine the nb of spawned tasks (xpvm and pvm> ok)
    - uses group for token ring communication
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include "pvm3.h"

#define GRPNAME "tokenring"

void matrix_load_par ( char nom[], double *tab, int N, int me, int nproc ) {
    FILE *f;
    int i,j,i_local;
    double *buffer = (double*)malloc(N*sizeof(double));
    i_local = 0;

    if(me==0){
      if ((f = fopen (nom, "r")) == NULL) { perror ("matrix_load_par : fopen "); } 
      for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fscanf (f, "%lf", (buffer+j));
            }
        if(i%nproc==0){
            memcpy(tab+i_local*N, buffer, N*sizeof(double));
            i_local++;
        }
        else{
            pvm_initsend( PvmDataDefault );
            pvm_pkdouble( buffer, N, 1 );
            pvm_send( pvm_gettid(GRPNAME, i%nproc), 0);
        }
      }
      fclose (f);
    }
    else{
        for(i=0; i<N/nproc; i++){
            pvm_recv( pvm_gettid(GRPNAME, 0), -1);
            pvm_upkdouble( tab+i*N, N, 1 );
        }
    }
}

void matrix_save_par ( char nom[], double *tab, int N, int nproc, int me ) {
    FILE *f;
    int i,j;
    double *buffer = (double*)malloc(N*sizeof(double));
    int i_local =0;

    if(me==0){
      if ((f = fopen (nom, "w")) == NULL) { perror ("matrix_save_par : fopen "); } 
      for (i=0; i<N; i++) {
        if(i%nproc==0){
            for (j=0; j<N; j++) {
                fprintf(f, "%8.2f ", *(tab+i_local*N+j));
            }
            i_local++;
        }
        else{
            pvm_recv( pvm_gettid(GRPNAME, i%nproc), -1);
            pvm_upkdouble(buffer, N, 1 );
            for(j=0; j<N; j++){
                fprintf(f, "%8.2f ", *(buffer+j)); 
            }
        }

        fprintf (f, "\n");
      }
      fclose (f);

    }else{
        for(i=0; i<N/nproc; i++){
            pvm_initsend( PvmDataDefault );
            pvm_pkdouble( tab+i*N, N, 1 );
            pvm_send( pvm_gettid(GRPNAME, 0), 0);
        }
        

    }
}

void matrix_display ( double *tab, int  N, int nproc ) {
  int i,j;

  for (i=0; i<N/nproc; i++) {
    for (j=0; j<N; j++) {
      printf ("%8.2f ", *(tab+i*N/nproc+j) );
    }
    printf ("\n");
  }

}

void gauss ( double * tab, int N ) {
  int i,j,k;
  double pivot;

  for ( k=0; k<N-1; k++ ){ /* mise a 0 de la col. k */
    /* printf (". "); */
    if ( fabs(*(tab+k+k*N)) <= 1.0e-11 ) {
      printf ("ATTENTION: pivot %d presque nul: %g\n", k, *(tab+k+k*N) );
      exit (-1);
    }
    for ( i=k+1; i<N; i++ ){ /* update lines (k+1) to (n-1) */
      pivot = - *(tab+k+i*N) / *(tab+k+k*N);
      for ( j=k; j<N; j++ ){ /* update elts (k) - (N-1) of line i */
	*(tab+j+i*N) = *(tab+j+i*N) + pivot * *(tab+j+k*N);
      }
      /* *(tab+k+i*N) = 0.0; */
    }
  }
  printf ("\n");
}

void gauss_par(double *tab, int N, int nproc, int me){
    int i, j, k, i_local, l, o;
    double pivot, m, n;


    i_local = 0;

    for(k=0;k<N-1;k++){
        
        if(me==k%nproc){
            i_local++;
            pivot = *(tab+((k-me)/nproc)*N+k);
            for(l=0; l<nproc; l++){
                pvm_initsend(PvmDataDefault);
                pvm_pkdouble(&pivot, 1, 1);
                pvm_send(pvm_gettid(GRPNAME, l), 0);
            }
        }
        else{
            pvm_recv(pvm_gettid(GRPNAME, k%nproc), -1);
            pvm_upkdouble(&pivot, 1, 1);
        }

        for(j=k+1;j<N;j++){
            if(me==k%nproc){
            m = *(tab+((k-me)/nproc)*N+j);
            for(l=0; l<nproc; l++){
                pvm_initsend(PvmDataDefault);
                pvm_pkdouble(&m, 1, 1);
                pvm_send(pvm_gettid(GRPNAME, l), 0);
            }
            }
            else{
                pvm_recv(pvm_gettid(GRPNAME, k%nproc), -1);
                pvm_upkdouble(&m, 1, 1);
            }
            for(i=i_local; i<N/nproc;i++){
                n = *(tab+i*N+k);
                *(tab+i*N+j) -= n/pivot*m;  //(tab+i*N+j) --> case courante
            } 
        }

        for(o=i_local; o<N/nproc; o++){
            *(tab+o*N+k) = 0;
        }

    }
   
}

void gauss_par_2(double *tab, int N, int nproc, int me){
    int i, j, k, i_local, l, o;
    double pivot, m, n;
    double* buffer = (double*)malloc(N*sizeof(double));


    i_local = 0;

    for(k=0;k<N-1;k++){
        
        if(me==k%nproc){

            i_local++;

            for(j=0; j<N;j++){
                buffer[j] = *(tab+((k-me)/nproc)*N+j);
            }

            //pivot = *(tab+((k-me)/nproc)*N+k);
            for(l=0; l<nproc; l++){
                pvm_initsend(PvmDataDefault);
                pvm_pkdouble(buffer, N, 1);
                pvm_send(pvm_gettid(GRPNAME, l), 0);
            }

            for(j=k+1; j<N; j++){
                for(i=i_local; i<N/nproc; i++){
                    n = *(tab+i*N+k);
                    m = buffer[j];
                    pivot = buffer[k];

                     *(tab+i*N+j) -= n/pivot*m;
                }            
            }

        }
        else{
            pvm_recv(pvm_gettid(GRPNAME, k%nproc), -1);
            pvm_upkdouble(buffer, N, 1);
        

        for(j=k+1;j<N;j++){
            for(i=i_local; i<N/nproc;i++){
                n = *(tab+i*N+k);
                m = buffer[j];
                pivot = buffer[k];
                *(tab+i*N+j) -= n/pivot*m;  //(tab+i*N+j) --> case courante
            } 
        }

        

    }
for(o=i_local; o<N/nproc; o++){
            *(tab+o*N+k) = 0;
        }

    }
   
}

main(int argc, char ** argv) {
    int NPROC = 8;		/* default nb of proc */
    int mytid;                  /* my task id */
    int *tids;                  /* array of task id */
    int me;                     /* my process number */
    int i;

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
/*           all the tasks are equivalent at that point                     */

     dowork(argc, argv, me, tids, NPROC );

     pvm_lvgroup( GRPNAME );
     pvm_exit();
}

/* Simple example passes a token around a ring */

dowork(int argc, char ** argv, int me, int tids[], int nproc ) {
/*     int token;
     int src, dest;
     int count  = 1;
     int stride = 1;
     int msgtag = 4;

     // Determine neighbors in the ring 
     src = pvm_gettid (GRPNAME, (me - 1 + nproc) % nproc );
     dest= pvm_gettid (GRPNAME, (me + 1) % nproc );
     if( me == 0 ) { 
        token = dest;
        pvm_initsend( PvmDataDefault );
        pvm_pkint( &token, count, stride );
        pvm_send( dest, msgtag );
        pvm_recv( src, msgtag );
     }
     else {
        pvm_recv( src, msgtag );
        pvm_upkint( &token, count, stride
 );
        pvm_initsend( PvmDataDefault );
        pvm_pkint( &token, count, stride );
        pvm_send( dest, msgtag );
     }
     printf("P%d (%x) token ring done\n", me, pvm_mytid());*/


/*-----------------------------------------------------------------------------*/

    int N, i, j, k;
    double *tab, pivot;
    char nom[255];
    FILE *f;
    struct timeval tv1, tv2;	/* for timing */
    int duree;
      
    if (argc != 3){
      printf ("Usage: %s <matrix size> <matrix name>\n", argv[0]);
      exit (-1);
    } 
    N = atoi ( argv[1] );
    strcpy(nom, argv[2]);
    if ( (tab=malloc(N*N/nproc*sizeof(double))) == NULL ) {
      printf ("Cant malloc %d bytes\n", N*N*sizeof(double));
      exit (-1);
    }
    gettimeofday( &tv1, (struct timezone*)0 );
        matrix_load_par ( nom , tab, N, me, nproc);
    gettimeofday( &tv2, (struct timezone*)0 );
    duree = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
    //printf ("loading time: %10.8f sec.\n", duree/1000000.0 );

    gettimeofday( &tv1, (struct timezone*)0 );
        gauss_par_2 ( tab, N , nproc, me);
    gettimeofday( &tv2, (struct timezone*)0 );
    duree = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
    printf ("computation time: %10.8f sec.\n", duree/1000000.0 );

    //sprintf ( nom+strlen(nom), ".result" );
    //matrix_save_par ( nom, tab, N, nproc, me);
    //matrix_display ( tab , N, nproc );
    }
