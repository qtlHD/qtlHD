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
  auto sexcols = phenames[0..3];
  auto dircols = phenames[3..$];

  auto addcovar = new double[][](pheno.length, pheno[0].length);
  foreach(i, row; pheno) {
    foreach(j, el; row) {
      addcovar[i][j] = el.value;
    }
  }

  auto bc = form_cross("BC");
  auto bc_sex = new bool[][](sexcols.length, pheno.length);
  foreach(i, sexcol; sexcols) {
    auto sexdir = get_sex_and_cross(bc, sexcol, "", phenames, pheno);
    bc_sex[i] = sexdir[0];

    if(sexcol=="bothsex") {
      assert(!all_female(bc_sex[i]));
      assert(!all_male(bc_sex[i]));
      assert(!all_same_sex(bc_sex[i]));
    }
    if(sexcol=="allfemale") {
      assert( all_female(bc_sex[i]));
      assert(!all_male(bc_sex[i]));
      assert( all_same_sex(bc_sex[i]));
    }
    if(sexcol=="allmale") {
      assert(!all_female(bc_sex[i]));
      assert( all_male(bc_sex[i]));
      assert( all_same_sex(bc_sex[i]));
    }

    auto covar1 = xchr_covar(bc, sexdir[0], sexdir[1]);
    auto covar2 = xchr_covar(bc, sexdir);
    auto covar3 = xchr_covar(bc, sexdir[0], sexdir[1]);
    auto covar4 = xchr_covar(bc, sexdir);
  }

  auto f2 = form_cross("F2");
  auto f2_sex = new bool[][][](sexcols.length, dircols.length, pheno.length);
  auto f2_dir = new int[][][][](sexcols.length, dircols.length, pheno.length, 1);
  foreach(i, sexcol; sexcols) {
    foreach(j, dircol; dircols) {
      auto sexdir = get_sex_and_cross(f2, sexcol, dircol, phenames, pheno);
      f2_sex[i][j] = sexdir[0];
      f2_dir[i][j] = sexdir[1];

      if(sexcol=="bothsex") {
        assert(!all_female(f2_sex[i][j]));
        assert(!all_male(f2_sex[i][j]));
        assert(!all_same_sex(f2_sex[i][j]));
      }
      if(sexcol=="allfemale") {
        assert(all_female(f2_sex[i][j]));
        assert(!all_male(f2_sex[i][j]));
        assert(all_same_sex(f2_sex[i][j]));
      }
      if(sexcol=="allmale") {
        assert(!all_female(f2_sex[i][j]));
        assert(all_male(f2_sex[i][j]));
        assert(all_same_sex(f2_sex[i][j]));
      }

      if(sexcol=="bothdir") {
        assert(!all_same_cross_direction(f2_dir[i][j]));
      }
      if(sexcol=="allforw") {
        assert(all_same_cross_direction(f2_dir[i][j]));
      }
      if(sexcol=="allback") {
        assert(all_same_cross_direction(f2_dir[i][j]));
      }

      auto covar1 = xchr_covar(f2, sexdir[0], sexdir[1]);
      auto covar2 = xchr_covar(f2, sexdir);
    }
  }
}