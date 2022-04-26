#!/usr/bin/perl

# 11 gev e- to fixed proton target, s= 21.5158
# Q2 = 2.5 GeV^2
# z = 0.4
# pt = 0.6
# x  0.05--0.65



my $Q2 = 2.5;

my $s = 21.5;
my $z = 0.4;
my $pt = 1.2;

for (my $j = 1; $j<=1;$j++){
    if ($j==1){
	$s = 21.5;
    }else{
	$s = 9000;
    }
    for (my $k=1;$k<=6;$k++){
#	$pt = $k * 0.2 - 0.1;
	if($k==6){
	    $Q2 = 7.;
	}else{
	    $Q2 = $k + 0.5;
	}
	for (my $l = 1; $l <=8; $l++){
	    $z = 0.3 + (2*$l-1) * 0.025;
	    open(outfile,">kine/kine_$j\_$k\_$l\.dat");
	    for ( my $i = 0; $i!=13;$i++){
		my $x = 0.05 + 0.05*$i;
		print outfile "$x $Q2 $z $pt $s\n";
	    }
	    close(outfile);
	    #system("./eic kine_$j\_$k\_$l\.dat neutron_$j\_$k\_$l\.dat");
	}

    }
}
