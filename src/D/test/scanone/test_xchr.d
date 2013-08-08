/**
 * test xchr
 */

module test.scanone.test_xchr;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.algorithm;
import std.math;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.marker;
import qtl.core.genotype;
import qtl.core.phenotype;
import qtl.core.scanone.xchr;
import qtl.core.hmm.cross;

unittest {
  writeln("Unit test " ~ __FILE__);

  auto phenames = ["bothsex", "allfemale", "allmale", "bothdir", "allforw", "allback"];
  auto phetable =[[   "0",       "0",         "1",       "0",       "0",        "1"  ],
                  [   "1",       "0",         "1",       "0",       "0",        "1"  ],
                  [   "0",       "0",         "1",       "1",       "0",        "1"  ],
                  [   "1",       "0",         "1",       "1",       "0",        "1"  ],
                  [   "0",       "0",         "1",       "0",       "0",        "1"  ],
                  [   "1",       "0",         "1",       "0",       "0",        "1"  ],
                  [   "0",       "0",         "1",       "1",       "0",        "1"  ],
                  [   "1",       "0",         "1",       "1",       "0",        "1"  ],
                  [   "0",       "0",         "1",       "0",       "0",        "1"  ],
                  [   "1",       "0",         "1",       "0",       "0",        "1"  ],
                  [   "0",       "0",         "1",       "1",       "0",        "1"  ],
                  [   "1",       "0",         "1",       "1",       "0",        "1"  ]];
  Phenotype pheno[][];
  foreach(line; phetable) {
    Phenotype[] ps = std.array.array(map!((a) { return set_phenotype(a); })(line));
    pheno ~= ps;
  }

  auto bc = form_cross("BC");
  auto bc_bothsex   = get_sex_and_cross(bc, "bothsex",   "", phenames, pheno);
  auto bc_allfemale = get_sex_and_cross(bc, "allfemale", "", phenames, pheno);
  auto bc_allmale   = get_sex_and_cross(bc, "allmale",   "", phenames, pheno);

  assert(!all_female(bc_bothsex[0]));
  assert( all_female(bc_allfemale[0]));
  assert(!all_female(bc_allmale[0]));

  assert(!all_male(bc_bothsex[0]));
  assert(!all_male(bc_allfemale[0]));
  assert( all_male(bc_allmale[0]));

  assert(!all_same_sex(bc_bothsex[0]));
  assert( all_same_sex(bc_allfemale[0]));
  assert( all_same_sex(bc_allmale[0]));

  auto f2 = form_cross("F2");
  auto sexcols = phenames[0..3];
  auto dircols = phenames[3..$];
  auto f2_sex = new bool[][][](sexcols.length, dircols.length, pheno.length);
  auto f2_dir = new int[][][][](sexcols.length, dircols.length, pheno.length, 1);
  foreach(i, sexcol; sexcols) {
    foreach(j, dircol; dircols) {
      auto sexdir = get_sex_and_cross(f2, sexcol, dircol, phenames, pheno);
      f2_sex[i][j] = sexdir[0];
      f2_dir[i][j] = sexdir[1];
    }
  }
}