#!/usr/bin/perl
#
# @File fixchains.pl.pl
# @Author gelpi
# @Created 02-dic-2019 19:22:32
#

use strict;

my $baseChId = 'A';
my $ter=0;
my $antResN=0;
my ($hisdata, $pdb) = @ARGV;
my $hid={};
open HIS, $hisdata;
while (<HIS>) {
	chomp;
	my ($lb, $num) = split;
	$hid->{$num} = 1;
}

open PDB, $pdb;
while (<PDB>) {
    next if /(SOL|NA|CL)/;
    if (/^TER/) {
        $baseChId='B';
    }
    my $atomN;
    my $elem;
    my $resN;
    if (/^ATOM/) {
        chomp;
        if (substr($_,12,1) ne ' ') {
            $atomN = substr($_,13,3).substr($_,12,1);
        } else {
            $atomN = substr($_,12,4);
        }
        if (/ILE/) {
            $atomN =~ s/ HD(.)/HD1$1/;
            $atomN =~ s/ CD / CD1/;
        }
        $atomN =~ s/ OC1/ O  /;
        $atomN =~ s/ OC2/ OXT/;
        $elem = substr($_,13,1);           
        $resN = substr($_,22,4)+0;
	if ($resN < $antResN) {
		$baseChId='B';
	}
	$antResN = $resN;
	if (/HIS/) {
		if ($hid->{$resN}) {
			$_ =~ s/HIS/HID/;
		} else {
			$_ =~ s/HIS/HIE/;
		}
	}
        print substr($_,0,12).$atomN.substr($_,16,5).$baseChId.substr($_,22,-1)."$elem\n";
    } else {
        print;
    }
}

#ILE CD CD1
#ILE HD1 HD11
#ILE HD2 HD12
#ILE HD3 HD13
#HIS HIE


