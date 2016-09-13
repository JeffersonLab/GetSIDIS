
setenv LHAPDF /work/halla/solid/yez/lhapdf
setenv LHAPATH ${LHAPDF}/share/LHAPDF
setenv CTEQPDF /work/halla/solid/yez/EIC_Gluon/CTEQ6/cteq-pdf-1.0.4

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${LHAPDF}/lib${CTEQPDF}/lib
setenv LIBRARY_PATH ${LIBRARY_PATH}:${LHAPDF}/lib${CTEQPDF}/lib
