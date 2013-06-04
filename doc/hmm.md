# Hidden Markov models (HMM)

## Introduction

A key issue in QTL analysis is the handling of missing genotype
information: either to consider positions between markers
(&ldquo;interval mapping&rdquo;) or to handle less-than-fully
informative markers.

We use a hidden Markov model (HMM) that connects the underlying
&ldquo;true&rdquo; genotypes to the observed marker genotypes.
The main task is to calculate probabilities for each possible true
genotype at each position along the genome, given the observed marker
genotypes.

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
genotype at a given marker.  `emit` is the only part that involves the
observed marker genotypes, and so is the only part that needs to deal
with the relationship between observed and true genotypes.

## Genotype encodings

When observed marker genotypes are read from a data file, they are
converted to a system for encoding
