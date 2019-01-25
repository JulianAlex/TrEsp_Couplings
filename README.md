# TrEsp_Coupling
  
  Julian Adolphs, 2017
 
 This Programs prepares the input-files for a calculation of Coulomb couplings 
 of pigments in pigment-protein-complexes by solving the Poisson-Boltzmann-Eqaution
 with the programm mead/multiflex (Bashford and Karplus). 
 
 Additionally the Coulomb couplings in vacuum are calculated. 
 
 ----------------------------------------------------------------------------
 
  "couplings" calculates:  

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

--------------------------------------------------------------------------------------------------

  Warning: If there are more than one pigment types in the input-file, Programm has to be adapted!  
           
--------------------------------------------------------------------------------------------------           
