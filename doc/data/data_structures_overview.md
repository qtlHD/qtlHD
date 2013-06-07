## Data structures for qtlHD

The data for an experimental cross for QTL analysis has a number of
parts.  Here we'll attempt to explain the basic pieces.

### Genotype data

For each individual or line, we will have genotypes at a set of
markers along chromosomes.

In the simplest cases (when there are two inbred founders), the
genotypes will form a two-dimensional array (i.e., a rectangle), with
rows corresponding to individuals (or lines) and columns corresponding
to markers.

In the more complex cases (e.g., Collaborative Cross, diversity
outcross, heterogeneous stock), we will have genotype data on the
individuals or lines and then another two-dimensional array with
genotypes for the founder lines.

### Marker maps

Typically, we will have a genetic map of the typed markers: an array
of chromosomes (including autosomes and possibly the X chromosome, and
note that the X chromosome may have an ID that is not `"X"`), with
each chromosome being an array of marker positions, in centiMorgans
(cM).

It will be important, in most cases, to also keep track of a physical
map of the markers, with positions in megabasepairs (Mb or Mbp).

While QTL analysis will be performed using the genetic map, it is
often important to translate the results to physical positions.
Similarly, we will want to translate the physical positions of things
like genes to genetic map positions.

### Phenotypes

In the simplest cases (e.g., backcross, intercross), the phenotypes
will correspond to the individuals that were genotyped.  So if the
genotype data are *n* &times; *m*, the phenotype data would be *n* &times;
*p*, with *n* = number of individuals, *m* = number of markers, *p* = number
of phenotypes.

But for recombinant inbred lines and the collaborative cross, we will
want to be able to separate the line-level genotypes and the
individual-level phenotypes.  That is, the genotypes will be for *n*
lines, but there will be *N* phenotyped individuals, with *N* &#8811;
*n*.  We will need a mapping from the *N* individuals to the *n*
lines.

Phenotypes will generally be numeric, and mostly doubles, though they
could be ints (e.g., 0 or 1, for a binary trait), though ints might be
treated as doubles in the qtlHD software, for simplicity.

In some cases, though, it will be convenient to include character
strings as phenotypes. For example, this is useful for tracking the order of
crosses among founder lines to generate a Collaborate Cross line.  In
additional, character strings can be useful for describing different
environments or treatments.  We might just require that these be
encoded as integers, or we might want a separate table of non-numeric
phenotypes.

In R/qtl, we had included all phenotypes and covariates (e.g.,
environment, sex, and treatment) together.  They might be separated,
but I see little advantage to that.

There will be particular &ldquo;phenotypes&rdquo; of particular
importance for making sense of the genotype data: basically sex and
cross direction (or other individual-level information about cross,
such as the number of generations in an advanced intercross
population).

### Meta-covariates

It will be important to include another table of covariates describing
the phenotypes (or covariates), which I like to call
&ldquo;meta-covariates&rdquo;.  For example, for phenotypes that
correspond to the expression of genes, we would want gene symbols and
genomic positions, and perhaps tissue.  For a phenotype measured over
time, we would want the time that corresponds to each phenotype
column.

Suppose the phenotype data are *n* &times; *p*, with *n* = number of
individuals and *p* = number of phenotypes.  Shat is needed is another
table that is *p* &times; *q*, where *q* = number of meta-covariates.
Much of this table will be missing.  For example, if there is gene
expression data plus body size measured over
time, there might be six meta-covariates: an indicator for whether a
phenotype column corresponds to gene
expression data; an indicator for whether a phenotype column
corresponds to body size; chromosome, position, and gene symbol for
the gene expression phenotypes; and time for the body size phenotypes.
Chromosome, position, and gene symbol would be missing for the body size
phenotypes, while time might be missing for the gene expression
phenotypes.


