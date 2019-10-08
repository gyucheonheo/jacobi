#include <sys/time.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define ERROR -1
double getError(double **, double *, double *, int);

int flag=0;
void
usage(){
  printf("\n");
  printf("USAGE: ./jacobi_mp -t <threads> -n <size> -i <iteration> -c <convergence>\n");
  printf("\n");
  printf("OPTIONS: \n");
  printf("\t -t Integer \t\t The number of threads\n");
  printf("\t -n Integer \t\t The size of n x n matrix for A\n");
  printf("\t -i Integer \t\t The number of iterations \n");
  printf("\t -c Double  \t\t The threshold of convergence \n");
    
  printf("Example : \n");
  printf("   ./jacobi_mp -t 4 -n 5 -i 25 -c 0.00001\n\n");
}  

void
jacobi_free(double ** arr, int n){
    int i;
    for(i = 0; i < n; i++){
        free(arr[i]);
    }
    free(arr);
}

void
getDiagonal(double **A, double **D, int n){
    
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n ; i++){
        for(int j = 0; j < n ; j++){
            if (i == j){
                D[i][i] = A[i][i];
            } else{
                D[i][j] = 0.0;
            }
        }
    }
}

void
getRemainder(double **A, double **R, int n){
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if (i == j){
                R[i][j] = 0.0;
            } else {
                R[i][j] = A[i][j];
            }
        }
    }
}

static void
pprintf(double **array, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lf ", array[i][j]);
        }
        printf("\n");
    }
}

double*
jacobi(double **A, double *b, double *x, int iteration, int n, double convergence, int thread){
    double **R = malloc( n * sizeof(double *));
    double **D = malloc( n * sizeof(double *));
    double *x_temp;
    /* Every task is independent */
    #pragma omp parallel for
    for(int k = 0; k <  n; k++){
        R[k] = (double *)malloc(n * sizeof(double) );
        D[k] = (double *)malloc(n * sizeof(double) ); 
    }
    
    getDiagonal(A, D, n);
    getRemainder(A, R, n);
    
    /* Any way to make progress on this for-loop ? as it is in charge of "iteration" */
    
    double *prev_x;
    for(int i = 0 ; i < iteration && !flag ; i++){
    #pragma omp parallel for num_threads(thread) shared(flag)
        for(int p = 0; p < n; p++){
            double result = 0.0;
            for(int q = 0; q < n; q ++){
                result += (R[p][q]*x[q]);
            }
           x[p] = (b[p] - result)/D[p][p];
       }
       if (i == 0){
          prev_x = x;
       }
      double tolerance = 0.0;
      double tolerance_c = 0.0;
      double tolerance_p = 0.0;

      if ( i > 0 && i < n ){
          for(int xi = 0; xi < n; xi++){
                    tolerance_c += pow(x[xi]-prev_x[xi], 2.0);
                    tolerance_p += pow(x[xi], 2.0);
                 }
          tolerance = sqrt(tolerance_c) / sqrt(tolerance_p);
          if ( (i > iteration) && (tolerance < convergence) ){ 
                     flag = 1;
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
int main(int argc, char **argv){
   int opt;
   int thread;
   int tflag=0;
   int n;
   int nflag=0;
   int i;
   int iflag=0;
   double c;
   int cflag=0;
   int test = 0;
   while ( (opt = getopt(argc, argv, "t:n:i:c:d")) != -1){
       switch(opt){
        case 't':
          thread = atoi(optarg);
          tflag = 1;
         break;
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
       }
   }
  if ( (tflag & nflag & iflag & cflag) == 0 ){
    usage();
    exit(ERROR);
  }

  double **A = malloc(n * sizeof(double *));
  double *b = malloc( n * sizeof(double));
  double *x = malloc( n * sizeof(double));

  for(int k = 0; k < n; k++){
     A[k] = (double *)malloc(n*sizeof(double)); 
   }

  if ( test == 1){  
      A[0][0] = 2;
      A[0][1] = 1;
      A[1][0] = 5;
      A[1][1] = 7;

      b[0] = 11;
      b[1] = 13;

      x[0] = 1;
      x[1] = 1;
    }
  else {
      int max_row = n*n; 
      for(int k=0; k < n; k++){
          b[k] = rand()% max_row; 
          for(int l=0; l < n; l++){
            if (k==l){
                A[k][l] = rand()%(2*n) + max_row;
            } else{
                A[k][l] = rand() % n;
            }
        }
    }
  }


   
  double start = omp_get_wtime();
  x = jacobi(A, b, x, i, n, c, thread);
  double end = omp_get_wtime();
  double error = getError(A, x, b, n);
  printf("Elapsed time : %lf\n", end-start);
  printf("Error : %lf\n", error);
  jacobi_free(A, n);
  return 0;
}

double
getError(double **A, double *x, double *b, int n){
    double error;
    double val;
    for(int i = 0; i < n ; i++){
        val = 0.0;
        for(int j = 0; j < n; j++){
            val += A[i][j] * x[j];
        }
        val -= b[i];
        error += pow(val, 2.0);
    }
    return sqrt(error);
}
