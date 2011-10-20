#!/usr/bin/perl

open(OUT1, ">../map.csv") or die("Cannot write to map file\n");
open(OUT2, ">../parents.csv") or die("Cannot write to parents file\n");
print OUT1 "marker,chr,pos\n";
print OUT2 "marker,A/J,C57Bl/6J,129SvlmJ,NOD/LtJ,NZO/HlLtJ,CAST/EiJ,PWK/PhJ,WSB/EiJ\n";

foreach $i (18..19, "X") {
    if($i ne "X" and $i < 10) { 
	$ifile = "chr0$i\.csv"; 
    }
    else { 
	$ifile = "chr$i\.csv"; 
    }

    open(IN, $ifile) or die("Cannot read from $ifile"); 
    $line = <IN>; $line = <IN>;
    while($line = <IN>) {
	chomp($line);
	@v = split(/,/, $line);
	if($v[1] eq "20") { $v[1] = "X"; }

	@g = @v[6..13];
	foreach $g (@g) {
	    if($g ne "A" and $g ne "C" and $g ne "G" and $g ne "T") {
		print "odd allele: $v[0] $g\n";
	    }
	}

	print OUT1 (join(",", @v[0..2]), "\n");
	print OUT2 (join(",", @v[(0, 6..13)]), "\n");

    }
}

