#include "couplings.h"    

/****************************************************************************************** 

Structure of pqr-format

   Field_name Atom_number Atom_name Residue_name Chain_ID Residue_number X Y Z Charge Radius

1   Field_name
      A string which specifies the type of PQR entry and should either be ATOM or HETATM. 
2   Atom_number
      An integer which provides the atom index.
3   Atom_name
      A string which provides the atom name.
4   Residue_name
      A string which provides the residue name.
  5 Chain_ID (skipped here)                                       
      An OPTIONA string which provides the chain ID of the atom. 
6   Residue_number
      An integer which provides the residue index.
7   X Y Z
      3 floats which provide the atomic coordiantes.
8   Charge
      A float which provides the atomic charge (in electrons).
9   Radius
      A float which provides the atomic radius (in Å).

******************************************************************************************/ 



void outPqrFile(char *in_dat, int n1, int n3, int n4, char type[n1][5], int *n_atom, 
		  char atom_type[n1][4], char molec[n1][4], 
		  char chain[n1][2], int *n_res, double *x, double *y, double *z, double *q, 
		  double *r){  

  FILE *out_file;
  char *out_dat;

  char string[20], string1[7], string2[7], molec_name[7];
  char delim[] = ".";  // strtok cut after delimiter 
  int i, j, k;

  strcpy(string, in_dat); 
  strcpy(string2, ".pqr");
  out_dat = strtok(string, delim); // cuts file-ending, for example .pdb
  strcat(out_dat, string2);        // puts .pqr as ending
      
  out_file = fopen(out_dat, "w");
     
  printf("Achtung, Atomradien wurden skaliert mit Faktor %3.1lf !!! (couplings_output.c)\n\n", VDW_FAC);

  for(i=0; i<n3; i++){  // loop over all pigments

    strcpy(string1, molec[0]); 
    strcpy(molec_name, string1);
    k = i*n4;

    for(j=0; j<n4; j++) 
      fprintf(out_file,"%s  %5d %4s %3s   %3d    %8.3lf %7.3lf %7.3lf %10.6lf  %6.3lf\n", 
	      type[j+k], n_atom[j+k], atom_type[j+k], molec_name, n_res[j+k], 
	      x[j+k], y[j+k], z[j+k], q[j+k], VDW_FAC*r[j+k]); 
  }

  fclose(out_file); 
  
}


//====================================================================================


int outMeadStFiles(int n1, int n3, int n4, char molec[n1][4], char atom_type[n1][4], 
		    double *x, double *y, double *z, double *q ){  

  FILE *out_file;
  char *out_dat;
  int j; 
  double q_0 = 0.0;
  char string1[7], string2[7];


  strcpy(string1, molec[0]); 
  strcpy(string2, ".st");
  out_dat = strcat(string1, string2);        
  printf("%s %s\n", string1, string2);

  out_file = fopen(out_dat, "w");

  fprintf(out_file,"7.0\n");   
  for(j=0; j<n4; j++)      
    fprintf(out_file,"%s   %3s %9.6lf  %3.1lf\n", molec[0], atom_type[j], q[j], q_0 ); 
    
  fclose(out_file); 

  return 0; 

}

//==================================================================================

int outMeadSitesFiles(char *in_dat, int n1, int n3, char molec[n1][4]){  

  FILE *out_file;
  char *out_dat;

  char string[20], string2[7]; 
  char delim[] = ".";               // strtok cuts after delimiter 
  char string1[4]; 
  int i;

  strcpy(string, in_dat);           // pdb_file_name.pdb
  strcpy(string1, molec[0]);        // for example CLA
  strcpy(string2, ".sites");
  out_dat = strtok(string, delim);  // cuts file-ending, for example .pdb
  strcat(out_dat, string2);         // puts .sites as ending
      
  out_file = fopen(out_dat, "w");
     
  for(i=0; i<n3; i++){
    fprintf(out_file,"%d   %3s \n", i+1, string1 );   
  }
  fclose(out_file); 
  
  return 0; 

}

//==============================================================

int outGmFiles(char *in_dat){  

  FILE *out_file;
  char *out_dat;

  char string[20], string2[7];
  char delim[] = ".";                // strtok cuts after delimiter 

  strcpy(string, in_dat); 
  strcpy(string2, ".mgm");
  out_dat = strtok(string, delim);   // cuts file-ending, for example .pdb
  strcat(out_dat, string2);          // puts .mgm as ending
      
  out_file = fopen(out_dat, "w");
     
  fprintf(out_file,"ON_GEOM_CENT %d %6.3lf\n", GRIDSIZE1, GRIDRES1);     // edit in couplings.h !!
  fprintf(out_file,"ON_GEOM_CENT %d %6.3lf\n", GRIDSIZE2, GRIDRES2); 
  fprintf(out_file,"ON_GEOM_CENT %d %6.3lf\n", GRIDSIZE3, GRIDRES3); 
    
  fclose(out_file); 
  
 
  strcpy(string, in_dat); 
  strcpy(string2, ".ogm");
  out_dat = strtok(string, delim);   // cuts file-ending, for example .pdb
  strcat(out_dat, string2);          // puts .ogm as ending
      
  out_file = fopen(out_dat, "w");
     
  fprintf(out_file,"ON_GEOM_CENT %d %6.3lf\n", GRIDSIZE2, GRIDRES2); 
  fprintf(out_file,"ON_GEOM_CENT %d %6.3lf\n", GRIDSIZE3, GRIDRES3); 

  fclose(out_file); 
  

  return 0; 

}

//==============================================================





 
