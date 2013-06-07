## CSV genotype encodings

CSV files containing genotype, phenotype, and marker map information
can have varying codes for genotypes and for missing values.

The qtlHD function to read in data from a csv file should take as
arguments:

    file name
    cross type
    genotype codes
    missing value codes
    
These could each be character strings.  I would suggest that,
initially, we only allow four possible cross types: BC, F2, RISELF,
RISIB.  The genotype codes and missing value codes could be character
strings with white space separating the separate codes.

We could have defaults for the genotype and missing value codes,
specific for each cross type.

I would allow multiple possible missing value codes, separated by
white space.  I would take the default to be:

    missing_values = "NA -"

Though one might want `missing_values = "NA - n/a na N/A"`.

### BC

For a backcross, we need genotype codes for the homozygote and the
heterozygote.  The default would be:

    genotypes = "A H"

Another typical example would be `genotypes = "AA AB"`.

### F2

For an intercross, we need *five* genotype codes: for AA, AB, BB,
&ldquo;not BB&rdquo;, and &ldquo;not AA&rdquo;.  The default would be:

    genotypes = "A H B D C"

The order of these is a bit confusing, but this choices matches
R/qtl.

I suppose one possible problem here is that I'm assuming the
genotype codes do not contain spaces.  You could imagine someone using
`"not AA"` and `"not BB"`, but that wouldn't be allowed if the
genotype codes are provided as a character string with white space
between the codes.


### RISELF and RISIB

For recombinant inbred lines, we need two genotype codes, for the two
homozygotes.  The default would be:

    genotypes = "A B"
    
But another typical choice might be `genotypes = "AA BB"`.
