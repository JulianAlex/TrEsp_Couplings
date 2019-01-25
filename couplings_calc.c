#include "couplings.h"

void attachAtomProperties(int n1, int n2, char atom_type[n1][4], double *q, char atom_tc[n2][4], double *trans_char){       

  // attaches every atom of the pdb-file the corresponding charge/radii from the 
  // transition-charges/ atom-radii input-file
  // depending on atom-names NOT order or number!! 

  int i,j; 


  for(i=0; i<n1; i++){      // loop over all atoms in pdb-file
    for(j=0; j<n2; j++){     // loop over all atoms in transition-charge file
                                                     // compare atom names 
      if( strcmp( atom_type[i], atom_tc[j] ) == 0)   // if strings are equal => 0 
	q[i] = trans_char[j];
    }
  }


}


//======================================================================


void dipoleStrength(int n1, int n2, int n3, int n4, double *dipst, double *x, double *y, double *z, 
		    double **dpos, double **qx, double **qy, double *q, char atom_type[n1][4]){       

  // Partial_Charges from TrEsp-Method
  // Calculate dipolestrength for each chlorophyll 
  // Unit of Dipolstrength: e*Ang = 4.8 Debye

  double dx[n1], dy[n1], dz[n1];
  double abs_x, abs_y, abs_y_, beta, cosbeta;
  double xa, xb, xc, xd, ya, yb, yc, yd, za, zb, zc, zd;
  int    i, j, k, m;

  double *qy_1 = calloc( n3, sizeof(double) );
  double *qy_2 = calloc( n3, sizeof(double) );
  double *qy_3 = calloc( n3, sizeof(double) );



  for(i=0; i<n1; i++)
    q[i] *= RENORM;   // Renormalization of partial charges
		

  for(j=0; j<n3; j++){  // sum over all chlorophylls 

    k = j*n4;

    for(i=0; i<n4; i++){
      dx[i] = DEBYE*q[i+k]*x[i+k];
      dy[i] = DEBYE*q[i+k]*y[i+k];
      dz[i] = DEBYE*q[i+k]*z[i+k];
    }
    
    for(i=0; i<n4; i++){
      qy_1[j] += dx[i];
      qy_2[j] += dy[i];
      qy_3[j] += dz[i];
    }

    dipst[j] = sqrt(SQ(qy_1[j])+SQ(qy_2[j])+SQ(qy_3[j]));

    printf("%3d %10.6lf  dipole strength/Debey\n", j+1, dipst[j]);
  
  }

  printf("\n");

  //-----------------------------------------------------------------

  // Calculate Angle between dipole moment from partial charges and 
  // NB-ND axis (used for extended dipole and point dipole approx)

  
  char string1[4], string2[4], string3[4], string4[4];

  strcpy(string1,"NA"); strcpy(string2,"NB"); 
  strcpy(string3,"NC"); strcpy(string4,"ND");  
  


  for(i=0; i<n3; i++){     // loop over all pigments
    k = i*n4;

    for(j=0; j<n4; j++){   // loop over all atoms of pigment i

      if( strcmp( atom_type[j+k], string1 ) == 0){ 
	xa = x[j+k]; ya = y[j+k]; za = z[j+k]; 
      }
      if( strcmp( atom_type[j+k], string2 ) == 0){ 
	xb = x[j+k]; yb = y[j+k]; zb = z[j+k];
      }
      if( strcmp( atom_type[j+k], string3 ) == 0){ 
	xc = x[j+k]; yc = y[j+k]; zc = z[j+k];
      }
      if( strcmp( atom_type[j+k], string4 ) == 0){ 
	xd = x[j+k]; yd = y[j+k]; zd = z[j+k];
      }
      
    } // end j-loop


    qx[i][0]=xa-xc; qx[i][1]=ya-yc; qx[i][2]=za-zc;  // na-nc axis = qx-dipole of p.dip and ext.dip approx.	
    qy[i][0]=xb-xd; qy[i][1]=yb-yd; qy[i][2]=zb-zd;  // nb-nd axis = qy-dipole of p.dip and ext.dip approx.

    dpos[i][0]=(xb+xd)/2.0; dpos[i][1]=(yb+yd)/2.0; dpos[i][2]=(zb+zd)/2.0;  // qy-dipole position

    abs_x = sqrt(SQ(qx[i][0])+SQ(qx[i][1])+SQ(qx[i][2]));
    abs_y = sqrt(SQ(qy[i][0])+SQ(qy[i][1])+SQ(qy[i][2]));
      
    for(m=0; m<3; m++){
      qx[i][m] = qx[i][m]/abs_x;
      qy[i][m] = qy[i][m]/abs_y;
    }

    abs_y_ = sqrt(SQ(qy_1[i])+SQ(qy_2[i])+SQ(qy_3[i]));


    // Angle between NB-ND-axis and qy-direction from pointcharges

    cosbeta = (qy_1[i]*qy[i][0]+qy_2[i]*qy[i][1]+qy_3[i]*qy[i][2])/(abs_y_);  // qy[m] is normalized!
    beta = acos(cosbeta);
    
    printf("Angle between qy from transQ and NB-ND-axis:  "); 
    printf("cosbeta = %9.6f,  beta = %9.2f deg\n ", cosbeta, beta*180.0/PI); 
       
 
    // Angle between NA-NC-axis and qy-direction from pointcharges

    cosbeta = (qy_1[i]*qx[i][0]+qy_2[i]*qx[i][1]+qy_3[i]*qx[i][2])/(abs_y_);  // qx[m] is normalized!      
    beta = acos(cosbeta);
   
    printf("Angle between qy from transQ and NA-NC-axis:  ");
    printf("cosbeta = %9.6f,  beta = %9.2f deg\n ", cosbeta, beta*180.0/PI); 
    
  }// end i-loop


}


//==================================================================================================


void calcVacuumCouplings(int n3, int n4, double *x, double *y, double *z, double *q){       

  // Calculate Coulomb couplings with pointCharges on each atom
  // Partial_Charges from TrEsp-Method

  int i,j,k,m;
  double sum = 0.0, dx, dy, dz, dr;

  printf("\nCoulomb-couplings in vacuum for transition charges/cm-1:\n");
  printf("\n");

  for(i=0; i<n3; i++){      // sum over all chlorophylls 
    for(j=0; j<n3; j++){     // sum over all chlorophylls 
      if( j != i ){           // no self-interaction 
	sum = 0.0;
	for(k=0; k<n4; k++){    // sum over all atoms of a chlorophyll
	  for(m=0; m<n4; m++){   // sum over all atoms of a chlorophyll
	    dx = x[i*n4+k] - x[j*n4+m]; 
	    dy = y[i*n4+k] - y[j*n4+m];
	    dz = z[i*n4+k] - z[j*n4+m];
	    dr = sqrt( dx*dx + dy*dy + dz*dz ); // distance between to partial charges
	    sum += q[i*n4+k]*q[j*n4+m]/dr;

	    // printf("%d %d %d %d  %8.3lf  %8.3lf  %10.6lf  %10.6lf\n",  
	    //	   i, j, k, m, x[i*n4+k], x[j*n4+m], q[i*n4+k], q[j*n4+m]);

	  }
	} 
      }
      else if( j == i ){
	sum = 0.0;
      }
      sum *= 1000*ECHARGE*EVWZ/(4*PI*EPS_0*ANG); //cm-1
      printf("%3d %3d %10.3lf\n", i+1, j+1, sum);      
    }
    //printf("\n");
  }


}


//==================================================================================================

void dipDistance(int n3, double **dpos, double **qx, double **qy, 
		double **delta, double **rij){ 

  int i, j, k;
  double abs;

  // Calculate difference matrix between dipole-positions 
  // delta[i][j][k] => delta[i+j*nz][k] 

  for(i=0; i<n3-1; i++){ 
    for(j=i+1; j<n3; j++){ 
      for(k=0; k<3; k++)
	delta[i+j*n3][k] = dpos[i][k]-dpos[j][k];
    }
  }

  //  rij is abs.val. of distance (in Ang) between dipole_i and dipole_j  

    for(i=0; i<n3-1; i++){ 
      for(j=i+1; j<n3; j++){ 
	abs = 0;
	for(k=0; k<3; k++)
	  abs += SQ(delta[i+j*n3][k]);
	abs = sqrt(abs);
	rij[i][j] = abs;
	rij[j][i] = rij[i][j];
	for(k=0; k<3; k++)
	  delta[i+j*n3][k] = delta[i+j*n3][k]/abs;
      }
      rij[i][i]=0.0;
    }     
 
    for(i=0; i<n3-1; i++){ 
      for(j=i+1; j<n3; j++){ 
	for(k=0; k<3; k++)
	  delta[j+i*n3][k] = -delta[i+j*n3][k] ;
      }
    }

}

//==================================================================================================

void pdipApprox(int n3, double **qy, double **delta, double **vab, 
		 double **rij, double *dipst){

  // Coulomb interaction of two point-dipoles (approximation)  

  int i, j, k;
  double h0, h1, h2, hvek[n3][n3][3];


  for(i=0; i<n3-1; i++){ 
    for(j=i+1; j<n3; j++){ 
      h0=0; h1=0; h2=0;
      for(k=0; k<3; k++){
	h0 += delta[i+j*n3][k]*qy[i][k];
        h1 += delta[i+j*n3][k]*qy[j][k];
	h2 += qy[i][k]*qy[j][k];
      }
      hvek[i][j][0]=h0;
      hvek[i][j][1]=h1;
      hvek[i][j][2]=h2;
    }
  }	

  // Convert Energy [D^2/A^3] to E [cm-1],  Factor 5040 

  for(i=0; i<n3-1; i++){
    for(j=i+1; j<n3; j++){ 
      vab[i][j] = hvek[i][j][2]-3*hvek[i][j][1]*hvek[i][j][0];
      vab[i][j] *= DIELEK/(pow((rij[i][j]), 3))*5040.84*dipst[i]*dipst[j];
    }     
  }    
  
  for(i=0; i<n3; i++)
    vab[i][i] = 0.0;
  
  for(i=1; i<n3; i++){ 
    for(j=0; j<i; j++){ 
      vab[i][j] = vab[j][i]; 
    }
  }


  printf("\n\nCoulomb-couplings in vacuum for point dipole approximation/cm-1:\n");
  printf("\n");

  for(i=0; i<n3; i++){      // sum over all chlorophylls 
    for(j=i+1; j<n3; j++){     // sum over all chlorophylls 
      printf("%3d %3d %10.3lf\n", i+1, j+1, vab[i][j]);      
    }
  }
  printf("\n");

}


//============================================================================


void extDipApprox(int n3, double **qy, double **delta, double **vab, 
		  double **rij, double *dipst, double **dpos){

  // Coulomb interaction of two extended dipoles (approximation)  
  // Charges of extended dipole: q_i = dipst[i]/DEBYE/EXTENT,  unit[q_i]=e;  
  // dpos = center between NB-ND-coordinates, qy = direction of dipole

  double rmm,rmp,rpm,rpp;
  int i,j,k;

  for(i=0; i<n3; i++)
    for(j=0; j<n3; j++)
      vab[i][j] = 0.0;

  //  for(i=0; i<n3; i++)
  //    printf("Ext-Dip-Approx, l = %10.6lf Ang,  q = %5.3lf e,  mu=q*l= %5.3lf D\n",
  //	   EXTENT, dipst[i]/DEBYE/EXTENT, dipst[i]); printf("\n");

  for(i=0; i<n3-1; i++){ 
    for(j=i+1; j<n3; j++){ 
      rmm=0.0;rpm=0.0;rmp=0.0;rpp=0.0;
      for(k=0; k<3; k++){	
	rpp += SQ( (dpos[i][k]+0.5*EXTENT*qy[i][k]) - (dpos[j][k]+0.5*EXTENT*qy[j][k]) );
	rpm += SQ( (dpos[i][k]+0.5*EXTENT*qy[i][k]) - (dpos[j][k]-0.5*EXTENT*qy[j][k]) );
	rmp += SQ( (dpos[i][k]-0.5*EXTENT*qy[i][k]) - (dpos[j][k]+0.5*EXTENT*qy[j][k]) );
	rmm += SQ( (dpos[i][k]-0.5*EXTENT*qy[i][k]) - (dpos[j][k]-0.5*EXTENT*qy[j][k]) );
      }
      vab[i][j]  = (1.0/sqrt(rmm)-1.0/sqrt(rpm)-1.0/sqrt(rmp)+1.0/sqrt(rpp))*dipst[i]*dipst[j]/SQ(DEBYE*EXTENT) ;
      vab[i][j] *= 1000*DIELEK*ECHARGE*EVWZ/(4*PI*EPS_0*ANG);
    }     
  }    

  for(i=0; i<n3; i++)
    vab[i][i] = 0.0;
  
  for(i=1; i<n3; i++){ 
    for(j=0; j<i; j++){ 
      vab[i][j] = vab[j][i]; 
    }
  }

  printf("\nCoulomb-couplings in vacuum for extended dipole approximation/cm-1:\n");
  printf("\n");

  for(i=0; i<n3; i++){       // sum over all chlorophylls 
    for(j=i+1; j<n3; j++){     
      printf("%3d %3d %10.3lf\n", i+1, j+1, vab[i][j]);      
    }
  }
  printf("\n");

}





