#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define ERROR -1

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


int main(int argc, char **argv){
   int opt;
   int thread;
   int tflag=0;
   int n;
   int nflag=0;
   double c;
   int cflag=0;
   int i;
   int iflag=0;
   while ( (opt = getopt(argc, argv, "t:n:i:c:")) != -1){
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
          if( sscanf(optarg, "%lf", c) == 0){
            usage();
            exit(ERROR);
          }
          cflag = 1;
          break;
        default:
          usage();
          exit(ERROR);
       }
   }
  if ( (tflag & nflag & cflag & iflag) == 0 ){
    usage();
    exit(ERROR);
  } 

}
