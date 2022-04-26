#!/usr/bin/perl

my $prefix = "kine/";

open(infile,"filelist");
while(<infile>){
    $filename = $_;
    chomp($filename);
    $input_file = $prefix + $filename;
    $filename =~ m/kine([\w,\.]+)$/;
    $output_file = "proton/proton$1";
    if (-e "$output_file"){
    }else{
	system("./eic $input_file $output_file");
    }
}
close(infile);
