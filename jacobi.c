#include <omp.h>
#include <sys/time.h>
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
  printf("\t -c Double  \t\t The convergence \n"); 
  printf("Example : \n");
  printf("   ./jacobi -n 5 -i 25 -c 0.001 \n");
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
 * Args : (Original, Diagonal, Iteration, Remainder, size)
 */
double *
jacobi(double **A, double *b, double *x, int iteration, int n, double convergence){


  double **R = malloc( n * sizeof(double *) );
  double **D = malloc( n * sizeof(double *) );
  double *x_temp;
  //  double convergence;
  int k;
 
  for(k = 0; k < n; k++){
    R[k] = (double *)malloc(n * sizeof(double) );
    D[k] = (double *)malloc(n * sizeof(double) );
  }
  
  getDiagonal(A, D, n); 
  getRemainder(A, R, n);
  
  double *prev_x = malloc(n * sizeof(double));
  for(int i = 0; i < iteration; i++){
    /* Dot product of remainder R and guess x */
    for(int p = 0; p < n; p++){
        double result = 0.0;
        for(int q = 0; q < n; q++){
            result += (R[p][q]*x[q]);
        }
        x[p] = (b[p]- result)/D[p][p];
    }
    if ( i == 0 ){
        prev_x = x; 
    }
    double tolerance = 0.0;
    double tolerance_c = 0.0;
    double tolerance_p = 0.0;    
    if ( i > 0 && i < n ){
        for(int xi = 0; xi < n; xi++){
           tolerance_c += pow(x[xi]-prev_x[xi], 2.0);
           tolerance_p += pow(x[xi], 2.0);
        
        tolerance = sqrt(tolerance_c) / sqrt(tolerance_p); 
        if ( (i > iteration) && (tolerance < convergence) ){
            return x;
        }
    }
    if (i > 0){
        prev_x = x;
    }
  }
      
  jacobi_free(R,n);
  jacobi_free(D,n);
  return x;
}
}
int
main(int argc, char **argv){
  /* Declaration for arguments */

  int opt;
  int n;
  int nflag = 0;
  int i;
  int iflag = 0;
  double c;
  int cflag = 0;
  int test = 0;
  /* Arguments */

  while( (opt = getopt(argc, argv, "n:i:c:d")) != -1 ){
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
      c = atof(optarg);
      cflag = 1;
      break;
    case 'd':
      test = 1;
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


  double **A = malloc( n * sizeof(double *) );
  double *b = malloc ( n * sizeof(double));
  double *x = malloc ( n * sizeof(double));
  
  int k;
  /* Allocate memory */
  for(k = 0; k < n; k++){
    A[k] = (double *)malloc( n * sizeof(double) );
  }
  /* Put random numbers into A and B */
  
 /* for(k = 0; k < n; k++){
    for(l = 0; l < n; l++){
      A[k][l] = randfrom(10,10);
    }
    [k] = randfrom(10,10);
  }*/
  if (test == 1){
    A[0][0] = 2;
    A[0][1] = 1;
    A[1][0] = 5;
    A[1][1] = 7;

    b[0] = 11;
    b[1] = 13;
    
    x[0] = 1;
    x[1] = 1;
  } else {
    int max_row = n*n;
    for(int k = 0; k < n ; k++){
        b[k] = rand()% max_row;
        for(int l = 0; l < n; l++){
            if (k == l){
                A[k][l] = rand()%(2*n) + max_row;
            } else{
                A[k][l] = rand() % n;
            }
        }
    }
  }

  /* Run Jacobi */
  double start = omp_get_wtime();
  x = jacobi(A, b, x, i, n, c);
  double end = omp_get_wtime();
  /* Compute Error */
  double error = getError(A, x, b, n);
  printf("Elapsed time : %lf\n", end-start);
  printf("Error : %lf\n", error);
  free(b);
  free(x);
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
