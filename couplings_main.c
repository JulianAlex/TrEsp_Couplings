/*************************************************************************
  
  Julian Adolphs, 2017
 
  Programm calculates:  

  - Coulomb couplings in vacuum, using the transition charges of 
    Madjet, Abdurahman, Renger, J. Phys. Chem. B 110, 2006  
    (Input-file "transChargesChla.txt")

  - Coulomb couplings in point dipole approximation in vacuum

  - Coulomb couplings in extended dipole approximation

  - Programm prepares input-files for multiflex: .pqr, .mgm, .ogm, .sites 
	
------------------------------------------------------------------------------

  compile program:      make

  run programm with:   ./couplings pdb_input_file.pdb 

------------------------------------------------------------------------------

  after that, run MEAD / MULTIFLEX with:  ./rum_mead.sh

------------------------------------------------------------------------------

  convert output units of couplings with: ./umrechnen mead_output.g

------------------------------------------------------------------------------

  Warning: If there are pigment types different from Chla in the input-file, 
           Programm has to be adapted!!  

*****************************************************************************/
                  
#include "couplings.h"      // contains macros & definitions   
 
       
int main(int argc, char *argv[ ]){         
  

  char *cmd, *in_dat_1, *in_dat_2, *in_dat_3; 
  int n1, n2, n3, n4, n5, nrow, hn;
  

  if(argc<2)
  {
    printf( "Missing argument\n" );
    exit(8);
  }
  if((cmd = malloc(strlen(argv[1])+1)) == NULL)
  {
    printf( "Out of memory\n" );
    exit(8);
  }
  in_dat_1 = argv[1];

  in_dat_2 = "trCharB3LYP_chl_a.txt";   
  in_dat_3 = "atomRadii.txt";         

  n1 = countInput(in_dat_1);  // n1 = number of atoms in pdb-file
  n2 = countInput(in_dat_2);  // n2 = number of atoms in transition-charge file
                              // n3 = number of chlorophylls
                              // n4 = number of atoms per chlorophyll
  n5 = countInput(in_dat_3);  // n5 = number of atoms in atom-radii file

  double *x = calloc( n1, sizeof(double));  // from file in_dat_1
  double *y = calloc( n1, sizeof(double));  // from file in_dat_1
  double *z = calloc( n1, sizeof(double));  // from file in_dat_1
  double *q = calloc( n1, sizeof(double));  // from file in_dat_2
  double *r = calloc( n1, sizeof(double));  // from file in_dat_3

  double *pdb_1  = calloc( n1, sizeof(double));  // don't know what these numbers are
  double *pdb_2  = calloc( n1, sizeof(double));  // so called them pdb 1 and 2

  int    *n_atom = calloc( n1, sizeof(int)); // from file in_dat_1
  int    *n_res  = calloc( n1, sizeof(int)); // from file in_dat_1

  char type[n1][5];       // ATOM, HETERO    // from file in_dat_1
  char atom_type[n1][4];  // C1B, C2B, ...   // from file in_dat_1
  char molec[n1][4];      // CLA, BCLA, ...  // from file in_dat_1
  char chain[n1][2];      // A, B, C, ...    // from file in_dat_1
  char molec_atom[n1][5]; // CLAC, CLAN, ... // from file in_dat_1


  char atom_tc[n2][4];           // C1B, C2B, ...    // from file in_dat_2
  char atom_rad[n5][4];          // C1B, C2B, ...    // from file in_dat_3
  double *trans_char = calloc( n2, sizeof(double));  // from file in_dat_2
  double *radii      = calloc( n5, sizeof(double));  // from file in_dat_3

  
  n3 = readInputPDB(in_dat_1, n1, type, n_atom, atom_type, molec, chain, n_res, x, y, z, pdb_1, pdb_2, molec_atom);
  n4 = n1/n3; 

  double *dipst = calloc( n3, sizeof(double) ); 

  double **qx   = calloc( n3, sizeof(double*) );
  double **qy   = calloc( n3, sizeof(double*) );
  double **dpos = calloc( n3, sizeof(double*) ); 
  double **rij  = calloc( n3, sizeof(double*) ); 
  double **vab  = calloc( n3, sizeof(double*) ); 
  nrow = n3;
  while(nrow--){                       
    qx[nrow]   = calloc(3, sizeof(double));
    qy[nrow]   = calloc(3, sizeof(double));
    dpos[nrow] = calloc(3, sizeof(double));
    rij[nrow]  = calloc(n3, sizeof(double));
    vab[nrow]  = calloc(n3, sizeof(double)); 
  }

  hn = SQ(n3);
  double **delta = calloc( hn, sizeof(double*));  
  nrow = hn;
  while(nrow--){
    delta[nrow] =    calloc(3, sizeof(double)); 
  }


  readAtomProperties(in_dat_2, n2, trans_char, atom_tc);  // read transition charges
  readAtomProperties(in_dat_3, n5, radii, atom_rad);      // read atom radii

  attachAtomProperties(n1, n2, atom_type, q, atom_tc, trans_char);
  attachAtomProperties(n1, n5, atom_type, r, atom_rad, radii);

  
  dipoleStrength(n1, n2, n3, n4, dipst, x, y, z, dpos, qx, qy, q, atom_type); 
  calcVacuumCouplings(n3, n4, x, y, z, q);
  dipDistance(n3, dpos, qx, qy, delta, rij); 
  pdipApprox(n3, qy, delta, vab, rij, dipst);
  extDipApprox(n3, qy, delta, vab, rij, dipst, dpos);


  outPqrFile(in_dat_1, n1, n3, n4, type, n_atom, atom_type, molec, chain, n_res, x, y, z, q, r); 
  outMeadStFiles(n1, n3, n4, molec, atom_type, x,y,z,q); 
  outMeadSitesFiles(in_dat_1, n1, n3, molec); 
  outGmFiles(in_dat_1); 


   
  free(x); x = NULL;
  free(y); x = NULL;
  free(z); z = NULL;
  free(q); q = NULL;
  free(r); r = NULL;
  free(pdb_1); pdb_1 = NULL;
  free(pdb_2); pdb_2 = NULL;
  free(n_res); n_res = NULL;
  free(radii); radii = NULL;
  free(dipst); dipst = NULL;
  free(n_atom); n_atom = NULL;
  free(trans_char); trans_char = NULL;

  nrow = hn;
  while(nrow--)
    free(delta[nrow]);
  free(delta); delta = NULL;

  nrow = n3; 
  while(nrow--){
    free(qx[nrow]);
    free(qy[nrow]);
    free(rij[nrow]);
    free(vab[nrow]);
    free(dpos[nrow]);
  }
  free(qx); qx = NULL;
  free(qy); qy = NULL;
  free(rij); rij = NULL;
  free(vab); vab = NULL;
  free(dpos); dpos = NULL;

  return 0;

}



