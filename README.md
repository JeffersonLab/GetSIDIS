
# SIDIS Cross Section Model 
## with nuclear-PDF implemented for SoLID and EIC
## -- Zhihong Ye, 07/26/2016

## How to use the generator:

* Install LHAPDF6, download the package and compile it by yourself. 
  Specify the path in Makefile or define it in your system PATH and LIBARAY

* Edit Makefile to include the correct path of compiled LHAPDF and ROOTSYS. Then compile the generator:
```js
   cd ./generator/
   make
```

* Change the configuration file if needed, e.g. "eic_c12_pion.dat", and run the code
```js
         ./GetSIDIS eic_c12_pion.dat, if the "FileNo" value in the file is not "0",
         ./GetSIDIS eic_c12_pion.dat N, for N=1, 2, ..., if the "FileNo" value in the file is "0"
```

## To Do
  *      Use a better sampling method when generating XS based on physical distributions
  *      Add absolute lumi

## Update Version in 07/11/2020:
 * Add EPPS16 and nNNPDFv2.0 nuclear PDF grids when calculating XS using Lsidis.h.
 Please download EPPS16 and nNNPDFv.20 grids from LHAPDF site and place them into the 
 LHAPDF/share/ folder. 
 * Make sure th grids for the selected nucleus target/beam is availalbe otherwise
 free proton-PDF set will be used without warning. Note that only a limited list of targets 
 is given in SIDIS.h::Init(). Please add more (and download the grids!) by editing this file.
 * To-Do: Define the list of nucleus targets/beams in a structure file for easily maintaining. 

## Update Version in 02/27/2018:
* Add Lsidis.h to calculate the SIDIS XS based on Tianbo Liu's method
* Used LHAPDF version 6.2 or above which you have to installed by yourself
* Remove any LHAPDF independent model, since Lsidis.h will do the job
* Remove (or won't include nCTEQ) since Lsidis.h will do the job
* Add two flags, isphy_hp and isphy_hm to determine if the events are in the good physical region


## Update Version in 09/16/2016:
* Add LHAPDF6 (replace LHAPDF5.8). Please install your won LHAPDF6. LHAPDF5.8 is still in the code but blocked-out.
* Keep the "Lite" version which only include EPS09 and CTEQ, and also the full version
* When choosing output LUND format files, generate PI+(K+) and PI-(K-) separately.
* When chossing output root files based on XS, save PI+(K+) and PI-(K-) into *_1.root and *_2.root separately
* Add a new configuraiton, name "SPECT", for spectrometer type configuration, e.g. HRS

## Update Version in 09/13/2016:
* Block the feature of calling "LHAPDF" model since we don't need it now. I will add the LHAPDF6 later.

## Update Version in 08/29/2016:
* Remove the separation of pip (kp) and pim(km) in the code since the generated
   file always contains both the info positive and negative particles. 
     e.g., in the future, only run "eic_c12_kaon.dat" to get both Pi+ and Pi- info
     and in the root files, their difference are just the cross sections and weight etc.

## Update Version in 08/15/2016
  *      All input parameters are defined in a input file, e.g. `eic_c12_pip.dat`.
         Please follow the format of the given examples to change the input values
         To run the run, simply type:
```js
         ./GetSIDIS eic_c12_pip.dat, if the "FileNo" value in the file is not "0",
         ./GetSIDIS eic_c12_pip.dat N, for N=1, 2, ..., if the "FileNo" value in the file is "0"
```
  *      Allow the generator to generate events uniformaly or based on the cross
                section shape. For the later case, a MAX cross section value is needed to be
         specified in the input file. I haven't figured out an automatic way to find
         the max cross sections.

## Current Version Note
 *      CTEQ-PDF is implemented; Please install cteq-pdf-1.0.4, a C++ libary to get CTEQ data
 *      EPS09 is also implemented; Please include "eps09_cxx" in the compiler
 *      LHAPDF is kept for comparison. You can choose this in "GetSIDIS.C" named "model". It if off by default
 *      In EIC collider mode, when calculation kinematic, I use the *average* mass and momentum,
             e.g., both quantities are divided by "A", so I can get correct x and Q2.
 *      In the Cross Section model, I contruct PDFs for different quarks (u, ubar, d, dbar, s, sbar),
             which can be either from EPS09, CTEQ or LHAPDf;
             then I obtained fragmentation functions for proton or neutrons;
             the "proton" and "neutron" SIDIS cross sections can be obtained by combining the PDF and fragmentation functions.
             Then for a nulceus (A,Z), the total differential cross section = Z*dxs_p + N*dxs_n
 *      The current EIC acceptance is a guess. Please change the following lines in GetSIDIS.C if you know better values:
       ```js
              if(config=="EIC" ){
                   Mom_Min_e = 0.5; Mom_Max_e = 3.*momentum_ele; //for Electrons
                   Mom_Min_h = 0.0; Mom_Max_h = 50.0;//for Pions

                   //Assuming a full 4PI acceptance
                   Th_Min_e = 0.0; Th_Max_e = 180.0;
                   Th_Min_h = 0.0; Th_Max_h = 180.0;
             }
       ```


# From the old version 
# SIDIS Cross Section Model 
## -- Zhihong Ye, 06/10/2014

## Note:
* This model is extracted fro Xin Qian's SIDIS generator "collider.C"
which is included as a reference.
* The XS model is coded in "SIDIS.h".
* A generator is given in "GetSIDIS.C"
* A simpler example of using the model is given in "example.C"
* LHAPDF is needed to be installed. Unpar "lhapdf.5.8.8.tar.gz",
and follow the instruction to install it. Specify the path in "SIDIS.h" or
define ${LHAPDF}=/yourpath/
* The output files include both the ROOT format files and LUND format files,
where, four root files will generated: 
```js
  "*_1_0.root", Q2=<10, pt<=1
  "*_2_0.root", Q2=<10, pt>1
  "*_3_0.root", Q2>10, pt<=1.
  "*_4_0.root", Q2>10, pt>1.
 ```
For SoLID config, chain the first two root files.
The LUND format file saves all events in one file.
* A older version of the model wrotten in FORTRAN is also given in "sidis_fortran/", as a reference.
* I held no responsiability to this code but you are welcome to discuss questions with me :->
