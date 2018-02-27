
Imin=0
Imax=400

for ((i=Imin; i<Imax; i++))
do
    #./gofarm eic_c12_pion.dat $i
    #./gofarm eic_d2_pion.dat $i
    #./gofarm eic_d2_pion_free.dat $i
    ./gofarm eic_c12_pion_free.dat $i
done
