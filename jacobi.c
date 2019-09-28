#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ERROR -1
double getError(double **, double *, double *, int);

void
usage(){
  printf("\n");
  printf("USAGE: ./jacobi -n <size> -i <iteration> -c <convergence>\n");
  printf("\n");
  printf("OPTIONS: \n");
  printf("\t -n Integer \t\t The size of n x n matrix for A\n");
  printf("\t -i Integer \t\t The number of iterations \n");
  printf("\t -c Double  \t\t The threshold of convergence \n");
    
  printf("Example : \n");
  printf("   ./jacobi -n 5 -i 25 -c 0.00001\n");
}

void
jacobi_free(double ** arr, int n){
  int i;
  for( i = 0; i < n; i++){
    free(arr[i]);
  }
  free(arr);
}

double
randfrom(double min, double max){
  double range = (max - min);
  double div = RAND_MAX / range;
  return (min + (rand() / div));
}

/* 
 * Name : getDiagonal 
 * Args : (from, to, size)
 */
void
getDiagonal(double ** A, double **D, int n){
  int i;
  int j;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if (i == j){
	D[i][i] = A[i][i];
      } else {
      D[i][j] = 0.0;
      }
    }
  }
}

/* 
 * Name : getRemainder
 * Args : (from, to, size)
 */
void
getRemainder(double ** A, double **R, int n){
  int i;
  int j;
  /* Initialize and Assign*/
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if ( i == j){
	R[i][j] = 0.0;
      } else {
	R[i][j] = A[i][j];
      }
    }
  }
}

/*
 * Name : jacobi
 * Args : (Original, Diagonal, Remainder, size)
 */
int
jacobi(double **A, double *b, double *x, double *x_guess, int iteration, double c, int n){


  double **R = malloc( n * sizeof(double *) );
  double **D = malloc( n * sizeof(double *) );
  double *temp;

  int i,j;

  double guess;
  double convergence;
  double tdiff;
  double diff = 0.0;
  int k;
  for(k = 0; k < n; k++){
    R[k] = (double *)malloc(n * sizeof(double) );
    D[k] = (double *)malloc(n * sizeof(double) );
  }

  getDiagonal(A, D, n); 
  getRemainder(A, R, n);


  /* Looping */
  k=0;
  do{
    for(i=0; i < n; i++){ /* row */
      guess = 0.0;
      for(j=0; j < n; j++){ /* column */
	if ( i != j ){
	  guess += (A[i][j]*x[j]);	  
	}
      }
      x_guess[i] = (b[i] - guess)/A[i][i];
    }
    /* Swapping */
    temp = x;
    x = x_guess;
    x_guess = temp;

    int t;
    for(t = 0; t < n; t++){
      tdiff = x_guess[t] - x[t];
      diff += pow(tdiff, 2.0);
    }
    k++;
  }
  while ( (k < iteration) && (sqrt(diff) > c ));

  for(k = 0; k < n; k++){
    free(R[k]);
    free(D[k]);
  }
  free(R);
  free(D);
}


int
main(int argc, char **argv){
  /* Declaration for arguments */

  int opt;
  int n;
  int nflag = 0;
  double c;
  int cflag = 0;
  int i;
  int iflag = 0;


  /* Arguments */

  while( (opt = getopt(argc, argv, "n:i:c:")) != -1 ){
    switch(opt){
    case 'i':
      i = atoi(optarg);
      iflag = 1;
      break;
    case 'n':
      n = atoi(optarg);
      nflag = 1;
      break;
    case 'c':
      if( sscanf(optarg, "%lf", &c) == 0){
	usage();
	exit(ERROR);
      };
      cflag = 1;
      break;
    default:
      usage();
      exit(ERROR);
      break;
    }
  }

  if ( (iflag & nflag & cflag) == 0){
    usage();
    exit(ERROR);
  }

  /* Allocate memory */


  double **A = malloc( 3 * sizeof(double *) );

  double *b = malloc ( n * sizeof(double));
  double *x = malloc ( n * sizeof(double));
  double *x_guess = malloc( n * sizeof(double));

  int k;
  int l;
  double rsum, val;
  for(k = 0; k < n; k++){
    A[k] = (double *)malloc(3 * sizeof(double) );
  }
  /* Put random number into A */
  srand(0);
  for(k = 0; k < n; k++){
    rsum = 0.0;
    x[k] = 0.0;
    for(l = 0; l < n; l++){
      val = rand() / (double)RAND_MAX;
      A[k][l] = val;
      rsum+=val;
    }
    A[k][k] +=rsum;
    b[k] = rand()/(double)RAND_MAX;
  }
  


  /* Run Jacobi */
  jacobi(A, b, x, x_guess, i, c, n);

  /* Compute Error */
  double error = getError(A, x, b, n);

  printf("Error value : %lf\n", error);
  
  free(b);
  free(x);
  free(x_guess);
  jacobi_free(A,n);

  
}

double
getError(double **A, double *x, double *b, int n){
  double error;
  double val;
  for(int i = 0; i < n; i++){
    val=0.0;
    for(int j = 0; j < n; j++){
      val += A[i][j] * x[j];
    }
    val -= b[i];
    error += pow(val, 2.0);
  }
  
  return sqrt(error);
}
