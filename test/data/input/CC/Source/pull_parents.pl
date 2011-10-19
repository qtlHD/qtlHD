#!/usr/bin/perl
# pull parental strain genotypes
# treat hets as missing and fill in all missing data

srand(76122977);

@parents = ("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/LtJ", "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ");
foreach $parent (@parents) {
    $ucparent = uc $parent;
    $isparent{$parent} = 1;
}

$ifile = "MUGAGtypes.csv";
open(IN, $ifile) or die("Cannot read from $ifile");
$line = <IN>; chomp($line);
@head = split(/,/, $line);
foreach $i (0..(@head-1)) {
    if($isparent{$head[$i]}) {
	$col{$head[$i]} = $i;
	push(@parcol, $i);
    }
}

#foreach $parent (@parents) {
#    print "$parent : $col{$parent}\n";
#}
#print(join("|", @parcol), "\n");

$ofile1 = "../parents.csv";
$ofile2 = "../map.csv";
open(OUT1, ">$ofile1") or die("Cannot write to $ofile1");
open(OUT2, ">$ofile2") or die("Cannot write to $ofile2");
print OUT1 "marker";
foreach $parent (@parents) {
    print OUT1 ",$parent";
}
print OUT1 "\n";
print OUT2 "marker,chr,pos\n";

while($line = <IN>) {
    chomp($line);
    @v = split(/,/, $line);

    # fill in missing data
    %freq = (); $n=$nmis=0;
    foreach $i (3..(@v-1)) {
	if($v[$i] ne "N" and $v[$i] ne "H") { ($freq{$v[$i]})++; $n++; }
	else { $nmis++; }
    }

    if($n==0) { next; } # no data on this marker 

    if($nmis > 0) {
	@alle = keys %freq;
	$na = @alle;
	$prob0 = $freq{$alle[0]} / $n;
	foreach $i (3..(@v-1)) {
	    if($v[$i] eq "N" or $v[$i] eq "H") { 
		$r = rand();
		if($r <= $prob0) {
		    $v[$i] = $alle[0];
		}
		else {
		    $v[$i] = $alle[1];
		    if($na < 2) { print "error: $v[0]\n"; }
		}
	    }
	}
    }
    

    if($v[1] eq "Y" or $v[1] eq "M") { next; }
    print OUT2 (join(",", @v[0..2]), "\n");
    print OUT1 (join(",", @v[(0, @parcol)]), "\n");
}

