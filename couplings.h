#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "couplings_func.h"   

// Size FMO with 8 pigments: 52 Ang

#define GRIDSIZE1  47 // 55    // Grid Size, odd number!
#define GRIDSIZE2 101 //101
#define GRIDSIZE3 281
#define GRIDRES1  5.0 //5.0    // Grid Resolution
#define GRIDRES2  1.0 //1.0
#define GRIDRES3  0.2  

#define VDW_FAC 1.4 // Scaling factor for Van der Waals radii 

#define ECHARGE 1.60217733*pow(10,-19)
#define PI      3.141592654
#define EPS_0   8.854187818*pow(10,-12)
#define EVWZ    8.06554  
#define ANG     pow(10,-10) 
#define MILLI   0.001
#define SQ(x) ((x) * (x))
#define DEBYE 4.8
#define DIP_VAC sqrt(21.0) // sqrt(21.0) vacuum dipole strength for Chla 
#define EXTENT 5.703125    // Extension of entended transition-dipole approx in Ang
#define RENORM  0.704076   // Renorm-Factor calculated for planar BCla and Majet-charges and vac_dipst 6.09 D

                            //0.796375 // Renorm-Factor calculated for planar Chla and Majet-charges and vac_dipst

#define BETA   0.0  //-7//0.0// Winkel zwischen qy und NB-ND, nur PDIP !!!
#define DIELEK 1.0  // 0.8 
