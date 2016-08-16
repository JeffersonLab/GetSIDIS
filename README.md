
===       SIDIS Cross Section Model
==       with nuclear-PDF implemented
==       for SoLID and EIC
==       -- Zhihong Ye, 07/26/2016

== To Do: ==
       1) Find the maximum XS value automatically. Now it has to be defined in the input file
       2) Update LHAPDF to the newest version. Now it is still 5.8
       3) Add nCTEQ15 as an alternative nPDF model.

====== Update Version in 08/10/2016: ======
       1) All input parameters are defined in a input file, e.g. input_c12_pip.dat.
         Please follow the format of the given examples to change the input values
         To run the run, simply type:
                  ./GetSIDIS input_c12_pip.dat, if the "FileNo" value in the file is not "0", or
                  ./GetSIDIS input_c12_pip.dat N, for N=1, 2, ..., if the "FileNo" value in the file is "0"
       2) Allow the generator to generate events uniformaly or based on the cross
         section shape. For the later case, a MAX cross section value is needed to be
         specified in the input file. I haven't figured out an automatic way to find
         the max cross sections.

====== Current Version Note: ========
       1), CTEQ-PDF is implemented; Please install cteq-pdf-1.0.4, a C++ libary to get CTEQ data
       2), EPS09 is also implemented; Please include "eps09_cxx" in the compiler
       3), LHAPDF is kept for comparison. You can choose this in "GetSIDIS.C" named "model". It if off by default
       4), In EIC collider mode, when calculation kinematic, I use the *average* mass and momentum,
             e.g., both quantities are divided by "A", so I can get correct x and Q2.
       5), In the Cross Section model, I contruct PDFs for different quarks (u, ubar, d, dbar, s, sbar),
             which can be either from EPS09, CTEQ or LHAPDf;
             then I obtained fragmentation functions for proton or neutrons;
             the "proton" and "neutron" SIDIS cross sections can be obtained by combining the PDF and fragmentation functions.
             Then for a nulceus (A,Z), the total differential cross section = Z*dxs_p + N*dxs_n
       6), The current EIC acceptance is a guess. Please change the following lines in GetSIDIS.C if you know better values:
              if(config=="EIC" ){
                   Mom_Min_e = 0.5; Mom_Max_e = 3.*momentum_ele; //for Electrons
                   Mom_Min_h = 0.0; Mom_Max_h = 10.0;//for Pions

                   //Assuming a full 4PI acceptance
                   Th_Min_e = 0.0; Th_Max_e = 180.0;
                   Th_Min_h = 0.0; Th_Max_h = 180.0;
             }


----------------From the old version ----------------
================================
= SIDIS Cross Section Model =
= -- Zhihong Ye, 06/10/2014 =
================================
Note:
1, This model is extracted fro Xin Qian's SIDIS generator "collider.C"
which is included as a reference.
2, The XS model is coded in "SIDIS.h".
1, A generator is given in "GetSIDIS.C"
2, A simpler example of using the model is given in "example.C"
3, LHAPDF is needed to be installed. Unpar "lhapdf.5.8.8.tar.gz",
and follow the instruction to install it. Specify the path in "SIDIS.h" or
define ${LHAPDF}=/yourpath/
4, The output files include both the ROOT format files and LUND format files,
where, four root files will generated:
"*_1_0.root", Q2=<10, pt<=1
"*_2_0.root", Q2=<10, pt>1
"*_3_0.root", Q2>10, pt<=1.
"*_4_0.root", Q2>10, pt>1.
For SoLID config, chain the first two root files.
The LUND format file saves all events in one file.
5, A older version of the model wrotten in FORTRAN is also given in "sidis_fortran/", as a reference.
6, I held no responsiability to this code but you are welcome to discuss questions with me :->
