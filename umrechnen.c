// Julian Adolphs
// convert multiflex-result (.g-file) from mead-units to meV and wavenumbers (cm-1)
// MolName.g contains site-site interactions in units of charge squared per length
//
// conversion factor  E[e/aA] -> E[meV]   1000*ECHARGE/(4*PI*EPS_0*ANG) = 14399.65173
// 
// compile:  gcc -lm umrechnen.c -o umrechnen
//
// execute:  ./umrechnen mead_output.g
// 
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <string.h>

#define ECHARGE 1.60217733*pow(10,-19)
#define PI      3.141592654
#define EPS_0   8.854187818*pow(10,-12)
#define EVWZ    8.06554                 // meV => cm-1
#define ANG     pow(10,-10)

#define mead2mev 14399.65173                       // mead-units => meV,   factor = 14399.65173 
#define mead2wz 1000*ECHARGE*EVWZ/(4*PI*EPS_0*ANG) // mead-units => cm-1,  factor = 116140.967014 


int main(int argc, char *argv[ ]){

  FILE *in_file;
  char *in_dat_1, *cmd;
  int n1=0, ch=0;

  if(argc<2)
  {
    printf( "Missing argument\n" );
    return 1;
  }

  if((cmd = malloc(strlen(argv[1])+1)) == NULL)
  {
    printf( "Out of memory\n" );
    return 1;
  }

  in_dat_1 = argv[1];
 
 // === Count input lines ===========================

  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);
    exit(8);
  }
  while (1) {            
    ch = fgetc(in_file);
    if (ch == '\n')
      ++n1;
    if (ch == EOF)
      break;
  }
  fclose(in_file);

  int i, a[n1], b[n1];
  double x[n1], y[n1];


  in_file = fopen(in_dat_1, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat_1);
    exit(8);
  }
  for(i=0; i<n1; i++){
    fscanf(in_file, "%d %d  %lf\n", &a[i], &b[i], &x[i]);   
  }
  fclose(in_file);
   
  printf("   i,   j,      E_MEAD,              E_meV,               E_cm-1\n");

  for(i=0; i<n1; i++)
    printf("%4d %4d   %18.6e %18.6lf %18.6lf\n", 
	   a[i], b[i], x[i], x[i]*mead2mev, x[i]*mead2wz);


}
