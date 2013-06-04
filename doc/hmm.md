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

I would recommend...
