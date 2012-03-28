/**
 * hmm_general: hmm functions for all cross types
 **/

module qtl.core.hmm.hmm_general;

import qtl.core.primitives;
import qtl.core.genotype;

// return vector of all possible true genotypes
TrueGenotype[] allTrueGeno(CrossType ct)
{
  switch(ct) {
  case BC:
    return [ new TrueGenotype("1,1"), new TrueGenotype("2,1") ];
  case F2:
    return [ new TrueGenotype("1,1"), new TrueGenotype("1,2"), new TrueGenotype("2,2") ];
  case RILself: case RILsib:
    return [ new TrueGenotype("1,1"), new TrueGenotype("2,2") ];
  }
}


// return vector of all possible true genotypes, phase-known case
new TrueGenotype[] allTrueGenoPK(CrossType ct)
{
  switch(ct) {
  case BC:
    return [ new TrueGenotype("1,1"), new TrueGenotype("2,1") ];
  case F2:
    return [ new TrueGenotype("1,1"), new TrueGenotype("1,2"), new TrueGenotype("2,1"), new TrueGenotype("2,2") ];
  case RILself: case RILsib:
    return [ new TrueGenotype("1,1"), new TrueGenotype("2,2") ];
  }
}

unittest {
  writeln("    unit test allTrueGeno for BC");
  assert(allTrueGeno(CrossType.BC) == [ new TrueGenotype("1,1"), new TrueGenotype("2,1") ]);
}

