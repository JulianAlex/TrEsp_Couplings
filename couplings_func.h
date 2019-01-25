int countInput(char *in_dat);

int readInputPDB(char *in_dat, int n1, char type[n1][5], int *n_atom, 
		  char atom_type[n1][4], char molec[n1][4], 
		  char chain[n1][2], int *n_res, double *x, double *y, double *z, double *pdb_1, 
		 double *pdb_2, char molec_atom[n1][5]);

void readAtomProperties(char *in_dat, int n2, double *trans_char, char atom_tc[n2][4]);

void attachAtomProperties(int n1, int n2, char atom_type[n1][4], double *q, char atom_tc[n2][4], double *trans_char);

void dipoleStrength(int n1, int n2, int n3, int n4, double *dipst, double *x, double *y, double *z, 
		    double **dpos, double **qx, double **qy, double *q, char atom_type[n1][4]); 

  void calcVacuumCouplings(int n3, int n4, double *x, double *y, double *z, double *q);

void dipDistance(int n3, double **dpos, double **qx, double **qy, 
		 double **delta, double **rij);

void pdipApprox(int n3, double **qy, double **delta, double **vab, 
		double **rij, double *dipst);

void extDipApprox(int n3, double **qy, double **delta, double **vab, 
		  double **rij, double *dipst, double **dpos);

void outPqrFile(char *in_dat, int n1, int n3, int n4, char type[n1][5], int *n_atom, 
		  char atom_type[n1][4], char molec[n1][4], 
		  char chain[n1][2], int *n_res, double *x, double *y, double *z, double *q, 
		double *r); 

int outMeadStFiles(int n1, int n3, int n4, char molec[n1][4], char atom_type[n1][4], 
		   double *x, double *y, double *z, double *q );
  
int outMeadSitesFiles(char *in_dat, int n1, int n3, char molec[n1][4]);  

int outGmFiles(char *in_dat);  
