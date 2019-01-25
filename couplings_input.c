#include "couplings.h"  

   
int countInput(char *in_dat){

  FILE *in_file;
  int ch, n=0;
 
  in_file = fopen(in_dat, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(8);
  }

  while (1) {            //zaehle zeilenzahl
    ch = fgetc(in_file);
    if (ch == '\n')
      ++n;
    if (ch == EOF)
      break;
  }
  
  fclose(in_file);

  printf("In_dat = %s,  n = %3d\n", in_dat, n);

  return n;
}

//===========================================================================


int readInputPDB(char *in_dat, int n1, char type[n1][5], int *n_atom, 
		  char atom_type[n1][4], char molec[n1][4], 
		  char chain[n1][2], int *n_res, double *x, double *y, double *z, double *pdb_1, 
		  double *pdb_2, char molec_atom[n1][5]){


  FILE *in_file;
  int i, n3;

  in_file = fopen(in_dat, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(8);
  }

  for(i=0; i<n1; i++){
    fscanf(in_file, " %s %d %s %s %s %d  %lf %lf %lf %lf  %lf %s\n", 
	   type[i], &n_atom[i], atom_type[i], molec[i], chain[i], &n_res[i], 
	   &x[i], &y[i], &z[i], &pdb_1[i], &pdb_2[i], molec_atom[i]); 
  }
  fclose(in_file);

  
  n3 = n_res[n1-1]; 
  
  printf("\n");
  printf("%s contains %d chlorophylls. Correct???\n", in_dat, n3);
  printf("A chlorophyll has %d atoms. Correct???\n", n1/n3);
  printf("\n");

  // In PDB-files there are no transition charges! 

  return(n3);

}


//===========================================================================


void readAtomProperties(char *in_dat, int n2, double *trans_char, char atom_tc[n2][4]){


  FILE *in_file;
  int i;

  in_file = fopen(in_dat, "r");
  if (in_file == NULL) {
    (void)printf("Can not open %s\n", in_dat);
    exit(8);
  }

  for(i=0; i<n2; i++){
    fscanf(in_file, " %lf %s\n", &trans_char[i], atom_tc[i]); 
  }
  fclose(in_file);

 
}





