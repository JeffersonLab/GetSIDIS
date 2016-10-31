
Imin=0
Imax=10

for ((i=Imin; i<Imax; i++))
do
    #./gofarm eic_c12_pion.dat $i
    #./gofarm eic_c12_kaon.dat $i
    #./gofarm eic_p_pion.dat $i
    #./gofarm eic_p_kaon.dat $i
    ./gofarm eic_d2_pion.dat $i
    #./gofarm eic_d2_kaon.dat $i
    #./gofarm hrs_h3_pi.dat $i
    #./gofarm hrs_he3_pi.dat $i
done
