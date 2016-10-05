Imin=0
Imax=100
for ((i=Imin; i<=Imax; i++))
do
    ./gofarm eic_c12_pion.dat $i
    #./gofarm eic_c12_kaon.dat $i
    #./gofarm eic_p_pion.dat $i
    #./gofarm eic_p_kaon.dat $i
done
