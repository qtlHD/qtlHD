# Hidden Markov models (HMM)

## Introduction

A key issue in QTL analysis is the handling of missing genotype
information: either to consider positions between markers
(&ldquo;interval mapping&rdquo;) or to handle less-than-fully
informative markers.

We use a hidden Markov model (HMM) to connect the underlying
&ldquo;true&rdquo; genotypes to the observed marker genotypes.
The main task is to calculate probabilities for each possible true
genotype at each position along the genome, given the observed marker
genotypes.

In the actual QTL analysis, we just deal with these probabilities for
the true genotypes; the observed marker data provide no further
information.

For a given type of cross, there is some fixed set of true genotypes,
and then another set of possible observed marker genotypes.  The
possible observed genotypes may be viewed as subsets of the true
genotypes.

For example, in an intercross, the possible true genotypes
are {AA, AB, BB}, and the possible observed genotypes are {&ndash;, A, H, B,
&ldquo;not BB&rdquo;, &ldquo;not AA&rdquo;}, with &ndash; indicating missing
data.  You can view each of these possible observed genotypes as a
subset of true genotypes:

    - = {AA, AB, BB}
    A = {AA}
    H = {AB}
    B = {BB}
    not BB = {AA, AB}
    not AA = {AB, BB}


The HMM has three components: `init` concerns the probabilities for
the true genotypes at any one position, `step` concerns the transition
probabilities for the true genotypes from one position to the next,
and `emit` concerns the relationship between the observed and true
genotype at a given marker.  `init` and `step` concern only with the true
genotypes.  `emit` is the only part that involves the
observed marker genotypes, and so is the only part that needs to deal
with the relationship between observed and true genotypes.

## Genotype encodings

When observed marker genotypes are read from a data file, they are
converted to a system for encoding the observed marker genotypes.

For crosses involving a pair of founder strains (backcross,
intercross, RIL, doubled haploids, haploids), the system for encoding
genotypes is simple and will be defined in the data file.

### The X chromosome

There are three tricky aspects to the treatment of the X chromosome in
QTL analysis.

1. One must take into account the sex of an individual and the
direction of the cross, e.g., in an intercross, whether the
individuals came from the cross (A &times; B) &times; (A &times; B) or
from the cross (B &times; A) &times; (B &times; A), where in each case
the mother (dam) is listed first and the father (sire) second.

2. One generally cannot distinguish hemizygous genotypes in males from
homozygous genotypes in females, and in the genotype data files the
hemizygous genotypes would often be encoded as homozygous.

3. In the actual QTL analysis with the X chromosome, one generally
needs to include some special covariates to handle cross direction and
sex.

So, for the X chromosome, the set of possible true genotypes is
different, and it can depend on the *sex* of an individual and also on
the cross direction that led to that individual.

And so also the possible observed marker genotypes depend on the sex
and possibly the cross direction for an individual.

#### Recommendation

For the encodings of both the true genotypes and the observed marker
genotypes, the X chromosome can be treated in *exactly* the same way
as the autosomes.  We just need some system for identifying the sex of
each individual as well as its cross information, and we need to
handle the observed marker genotype data somewhat specially.

I would suggest:

1. When we need the individuals' sexes, it should be included as a
phenotype, and we should have some simple way to identify which one.
For example, a parameter `sexname` in the data file indicates the name
of the phenotype that corresponds to sex, plus something like a
`Sex Symbol` table that indicates the symbols for `female` and
`male`, so that one might write `F`/`M`, `Female`/`Male`, `0`/`1`,
etc.

2. Similarly, we need a system for capturing each individuals' cross
direction information.  For maximal flexibility, I'm inclined to say
that this should be a character string.  It could be included in the
phenotypes, again with a parameter like `crossinfoname` in the data
file that indicates the name of the corresponding phenotype.

3. The genotype symbols for the X chromosome need to be allowed to be
different for each sex and each cross direction, and they are
different from those used on the autosomes.  For example, in an
intercross, `AA` might correspond to `0,0` for a female but `0,1` for
a male.

4. Alternatively, in a cross with *n* founders, we could allow let
genotype *n+1* correspond to the Y chromosome.


### Crosses with more than two founders

For crosses like the Collaborative Cross (CC), heterogeneous stock
(HS), the diversity outcross population (DO), and [MAGIC lines](http://www.ncbi.nlm.nih.gov/pubmed/19593375), in
which there are more than two founders, there is the added complexity
that the observed marker genotypes require two tables.  We would
generally have SNP genotypes, and we would need the SNP genotypes for
both the individuals phenotyped for QTL mapping and the SNP genoytypes
for the founder lines.

For each cross type, we would generally need individual-level
information on the cross direction, which might be encoded as a
character string.  Further, as with the simpler two-founder crosses,
we would need the sex of each individual.


## Not sure where to fit this in...

With *f* founders, there are *2<sup>f</sup>* possible phase-known
genotypes and *f* + *f(f-1)/2* possible phase-unknown genotypes.
When the X chromosome is involved, we have that many genotypes for
females and then another *f* genotypes for males.

We could fit the males into the current system for handling true
genotypes, but it would probably be best to expand it for males:  a
given genotype may be autosome or X-linked, and the X-linked ones can
be for a female (two alleles, just like autosome) or a male (just one
allele).

Alternatively, we could just let the *f* homozygous genotypes be
allowable for males, and use those to indicate that hemizygous
genotypes.

Yes, I think this final idea is the right one: identify the hemizygous
genotypes with the homozygous ones.
