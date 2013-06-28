/**
 * Unittests for cross.d
 **/

module test.hmm.cross_unittest;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

import qtl.core.primitives;
import qtl.core.genotype;
import qtl.core.hmm.cross;

unittest {
  writeln("Unit test " ~ __FILE__ ~ " BC");

  // unit test all_true_geno for BC, autosome
  auto bc = form_cross("BC");
  auto atg = bc.all_true_geno_A;
  assert(atg.length == 2);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,0));

  // blank stuff for autosome
  auto is_X_chr = false;
  auto is_female = false;
  auto cross_direction = new int[](1); // not used for BC

  // unit test init for BC
  foreach(gv; atg)
    assert(to!string(bc.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  // unit test step for BC
  double rf = 0.01;

  assert( to!string( bc.step(atg[0], atg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)) ));
  assert( to!string( bc.step(atg[0], atg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( bc.step(atg[1], atg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( bc.step(atg[1], atg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)) ));

  // unit test emit for BC
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeSymbolMapper("H");
  H ~= new TrueGenotype(1,0);

  double error_prob = 0.01;

  assert(bc.emit(NA, atg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(bc.emit(NA, atg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(to!string(bc.emit(A, atg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(A, atg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(bc.emit(H, atg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(H, atg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for BC
  assert(bc.nrec(atg[0], atg[0], is_X_chr, is_female, cross_direction) == 0.0);
  assert(bc.nrec(atg[0], atg[1], is_X_chr, is_female, cross_direction) == 1.0);

  assert(bc.nrec(atg[1], atg[0], is_X_chr, is_female, cross_direction) == 1.0);
  assert(bc.nrec(atg[1], atg[1], is_X_chr, is_female, cross_direction) == 0.0);
}


unittest {
  writeln("Unit test " ~ __FILE__ ~ " BC X chr female");

  // unit test all_true_geno for BC, X chr, female
  auto bc = form_cross("BC");
  auto atg = bc.all_true_geno_X;
  assert(atg.length == 3);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,0));
  assert(atg[2] == new TrueGenotype(1,1));

  // female X chr
  auto is_X_chr = true;
  auto is_female = true;
  auto cross_direction = new int[](1); // not used for BC

  // index to possible genotype
  auto ptgi = bc.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto ptg = new TrueGenotype[](ptgi.length);
  foreach(i, ptgival; ptgi)
    ptg[i] = atg[ptgival];
  assert(ptgi.length == 2);
  assert(ptgi[0] == 0);
  assert(ptgi[1] == 1);
  assert(ptg.length == 2);
  assert(ptg[0] == new TrueGenotype(0,0));
  assert(ptg[1] == new TrueGenotype(1,0));

  // unit test init for BC X chr female
  foreach(gv; ptg)
    assert(to!string(bc.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  // unit test step for BC X chr female
  double rf = 0.01;

  assert( to!string( bc.step(ptg[0], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)) ));
  assert( to!string( bc.step(ptg[0], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( bc.step(ptg[1], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( bc.step(ptg[1], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)) ));

  // unit test emit for BC X chr female
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeSymbolMapper("H");
  H ~= new TrueGenotype(1,0);

  double error_prob = 0.01;

  assert(bc.emit(NA, ptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(bc.emit(NA, ptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(to!string(bc.emit(A, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(A, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(bc.emit(H, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(H, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for BC X chr female
  assert(bc.nrec(ptg[0], ptg[0], is_X_chr, is_female, cross_direction) == 0.0);
  assert(bc.nrec(ptg[0], ptg[1], is_X_chr, is_female, cross_direction) == 1.0);

  assert(bc.nrec(ptg[1], ptg[0], is_X_chr, is_female, cross_direction) == 1.0);
  assert(bc.nrec(ptg[1], ptg[1], is_X_chr, is_female, cross_direction) == 0.0);
}


unittest {
  writeln("Unit test " ~ __FILE__ ~ " BC X chr male");

  // unit test all_true_geno for BC, X chr, male
  auto bc = form_cross("BC");
  auto atg = bc.all_true_geno_X;
  assert(atg.length == 3);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,0));
  assert(atg[2] == new TrueGenotype(1,1));

  // male X chr
  auto is_X_chr = true;
  auto is_female = false;
  auto cross_direction = new int[](1); // not used for BC

  // index to possible genotype
  auto ptgi = bc.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto ptg = new TrueGenotype[](ptgi.length);
  foreach(i, ptgival; ptgi)
    ptg[i] = atg[ptgival];
  assert(ptgi.length == 2);
  assert(ptgi[0] == 0);
  assert(ptgi[1] == 2);
  assert(ptg.length == 2);
  assert(ptg[0] == new TrueGenotype(0,0));
  assert(ptg[1] == new TrueGenotype(1,1));

  // unit test init for BC X chr male
  foreach(gv; ptg)
    assert(to!string(bc.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  // unit test step for BC X chr male
  double rf = 0.01;

  assert( to!string( bc.step(ptg[0], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)) ));
  assert( to!string( bc.step(ptg[0], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( bc.step(ptg[1], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( bc.step(ptg[1], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)) ));

  // unit test emit for BC X chr male
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto B = new GenotypeSymbolMapper("B");
  B ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(bc.emit(NA, ptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(bc.emit(NA, ptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(to!string(bc.emit(A, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(A, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(bc.emit(B, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(bc.emit(B, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for BC X chr male
  assert(bc.nrec(ptg[0], ptg[0], is_X_chr, is_female, cross_direction) == 0.0);
  assert(bc.nrec(ptg[0], ptg[1], is_X_chr, is_female, cross_direction) == 1.0);

  assert(bc.nrec(ptg[1], ptg[0], is_X_chr, is_female, cross_direction) == 1.0);
  assert(bc.nrec(ptg[1], ptg[1], is_X_chr, is_female, cross_direction) == 0.0);
}


unittest {
  writeln("Unit test " ~ __FILE__ ~ " F2");

  // unit test allTrueGeno for F2
  auto f2 = form_cross("F2");
  auto f2PK = form_cross_phaseknown(f2);
  auto g = f2.all_true_geno_A;
  auto gPK = f2PK.all_true_geno_A;

  assert(g.length == 3);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
  assert(g[2] == new TrueGenotype(1,1));

  // unit test allTrueGeno for F2 (phase-known)
  assert(gPK.length == 4);
  assert(gPK[0] == new TrueGenotype(0,0));
  assert(gPK[1] == new TrueGenotype(1,0));
  assert(gPK[2] == new TrueGenotype(0,1));
  assert(gPK[3] == new TrueGenotype(1,1));

  // blank stuff for autosome
  auto is_X_chr = false;
  auto is_female = false;
  auto cross_direction = new int[](1); // not used for autosome

  // unit test init for F2
  double val[] = [log(0.25), log(0.5), log(0.25)];

  foreach(i, gv; g)
    assert(to!string(f2.init(gv, is_X_chr, is_female, cross_direction)) == to!string(val[i]));

  // unit test init for F2 (phase-known)
  foreach(gv; gPK)
    assert(to!string(f2PK.init(gv, is_X_chr, is_female, cross_direction)) == to!string(val[0]));

  //  unit test step for F2
  double rf = 0.01;

  assert( to!string( f2.step(g[0], g[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2.step(g[0], g[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( f2.step(g[0], g[2], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( f2.step(g[1], g[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step(g[1], g[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf^^2 + (1-rf)^^2) ));
  assert( to!string( f2.step(g[1], g[2], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));

  assert( to!string( f2.step(g[2], g[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2.step(g[2], g[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( f2.step(g[2], g[2], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)^^2) ));


  // unit test step for F2 (phase-known)
  assert( to!string( f2PK.step(gPK[0], gPK[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[0], gPK[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[0], gPK[2], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[0], gPK[3], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( f2PK.step(gPK[1], gPK[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[1], gPK[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[1], gPK[2], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2PK.step(gPK[1], gPK[3], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( f2PK.step(gPK[2], gPK[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[2], gPK[2], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[2], gPK[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2PK.step(gPK[2], gPK[3], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( f2PK.step(gPK[3], gPK[3], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[3], gPK[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[3], gPK[2], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[3], gPK[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf^^2) ));

  // unit test emit for F2
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeSymbolMapper("H");
  H ~= new TrueGenotype(1,0);
  H ~= new TrueGenotype(0,1);
  auto B = new GenotypeSymbolMapper("B");
  B ~= new TrueGenotype(1,1);
  auto AorH = new GenotypeSymbolMapper("AorH");
  AorH ~= new TrueGenotype(0,0);
  AorH ~= new TrueGenotype(1,0);
  AorH ~= new TrueGenotype(0,1);
  auto HorB = new GenotypeSymbolMapper("HorB");
  HorB ~= new TrueGenotype(1,0);
  HorB ~= new TrueGenotype(0,1);
  HorB ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(f2.emit(NA, g[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2.emit(NA, g[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2.emit(NA, g[2], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2.emit(A, g[0], error_prob, is_X_chr, is_female, cross_direction))
         == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(A, g[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(A, g[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(H, g[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(H, g[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(H, g[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(B, g[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(B, g[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(B, g[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(AorH, g[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, g[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, g[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  assert(to!string(f2.emit(HorB, g[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(HorB, g[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(HorB, g[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));

  assert(f2PK.emit(NA, gPK[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2PK.emit(NA, gPK[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2PK.emit(NA, gPK[2], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2PK.emit(NA, gPK[3], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2PK.emit(A, gPK[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(A, gPK[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(A, gPK[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(A, gPK[3], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2PK.emit(H, gPK[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(H, gPK[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(H, gPK[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(H, gPK[3], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2PK.emit(B, gPK[3], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(B, gPK[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(B, gPK[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(B, gPK[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2PK.emit(AorH, gPK[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(AorH, gPK[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(AorH, gPK[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(AorH, gPK[3], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  assert(to!string(f2PK.emit(HorB, gPK[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2PK.emit(HorB, gPK[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(HorB, gPK[2], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(HorB, gPK[3], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob/2.0)));

  // unit test nrec for F2 (phase-known)
  assert( f2PK.nrec(gPK[0], gPK[0], is_X_chr, is_female, cross_direction) == 0.0 );
  assert( f2PK.nrec(gPK[0], gPK[1], is_X_chr, is_female, cross_direction) == 0.5 );
  assert( f2PK.nrec(gPK[0], gPK[2], is_X_chr, is_female, cross_direction) == 0.5 );
  assert( f2PK.nrec(gPK[0], gPK[3], is_X_chr, is_female, cross_direction) == 1.0 );

  assert( f2PK.nrec(gPK[1], gPK[0], is_X_chr, is_female, cross_direction) == 0.5 );
  assert( f2PK.nrec(gPK[1], gPK[1], is_X_chr, is_female, cross_direction) == 0.0 );
  assert( f2PK.nrec(gPK[1], gPK[2], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(gPK[1], gPK[3], is_X_chr, is_female, cross_direction) == 0.5 );

  assert( f2PK.nrec(gPK[2], gPK[0], is_X_chr, is_female, cross_direction) == 0.5 );
  assert( f2PK.nrec(gPK[2], gPK[2], is_X_chr, is_female, cross_direction) == 0.0 );
  assert( f2PK.nrec(gPK[2], gPK[1], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(gPK[2], gPK[3], is_X_chr, is_female, cross_direction) == 0.5 );

  assert( f2PK.nrec(gPK[3], gPK[3], is_X_chr, is_female, cross_direction) == 0.0 );
  assert( f2PK.nrec(gPK[3], gPK[1], is_X_chr, is_female, cross_direction) == 0.5 );
  assert( f2PK.nrec(gPK[3], gPK[2], is_X_chr, is_female, cross_direction) == 0.5 );
  assert( f2PK.nrec(gPK[3], gPK[0], is_X_chr, is_female, cross_direction) == 1.0 );
}


unittest {
  writeln("Unit test " ~ __FILE__ ~ " F2 X chr male");

  // unit test allTrueGeno for F2 X chr male
  auto f2 = form_cross("F2");
  auto f2PK = form_cross_phaseknown(f2);
  auto g = f2.all_true_geno_X;
  auto gPK = f2PK.all_true_geno_X;

  assert(g.length == 3);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
  assert(g[2] == new TrueGenotype(1,1));

  // unit test allTrueGeno for F2 X chr male (phase-known)
  assert(gPK.length == 4);
  assert(gPK[0] == new TrueGenotype(0,0));
  assert(gPK[1] == new TrueGenotype(1,0));
  assert(gPK[2] == new TrueGenotype(0,1));
  assert(gPK[3] == new TrueGenotype(1,1));

  // X chr male
  auto is_X_chr = true;
  auto is_female = false;
  auto cross_direction = new int[](1); // not used for male

  // index to possible genotype
  auto ptgi = f2.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto ptg = new TrueGenotype[](ptgi.length);
  foreach(i, ptgival; ptgi)
    ptg[i] = g[ptgival];
  assert(ptgi.length == 2);
  assert(ptgi[0] == 0);
  assert(ptgi[1] == 2);
  assert(ptg.length == 2);
  assert(ptg[0] == new TrueGenotype(0,0));
  assert(ptg[1] == new TrueGenotype(1,1));

  auto PKptgi = f2PK.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto PKptg = new TrueGenotype[](PKptgi.length);
  foreach(i, PKptgival; PKptgi)
    PKptg[i] = gPK[PKptgival];
  assert(PKptgi.length == 2);
  assert(PKptgi[0] == 0);
  assert(PKptgi[1] == 3);
  assert(PKptg.length == 2);
  assert(PKptg[0] == new TrueGenotype(0,0));
  assert(PKptg[1] == new TrueGenotype(1,1));

  // unit test init for F2 X chr male
  foreach(i, gv; ptg)
    assert(to!string(f2.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  // unit test init for F2 X chr male (phase-known)
  foreach(gv; PKptg)
    assert(to!string(f2PK.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  //  unit test step for F2 X chr male
  double rf = 0.01;

  assert( to!string( f2.step(ptg[0], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));
  assert( to!string( f2.step(ptg[0], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( f2.step(ptg[1], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( f2.step(ptg[1], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));

  // unit test step for F2 X chr male (phase-known)
  assert( to!string( f2PK.step(PKptg[0], PKptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));
  assert( to!string( f2PK.step(PKptg[0], PKptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( f2PK.step(PKptg[1], PKptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( f2PK.step(PKptg[1], PKptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));

  // unit test emit for F2 X chr male
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto B = new GenotypeSymbolMapper("B");
  B ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(f2.emit(NA, ptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2.emit(NA, ptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2.emit(A, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(A, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(B, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(B, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));

  assert(f2PK.emit(NA, PKptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2PK.emit(NA, PKptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2PK.emit(A, PKptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(A, PKptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  assert(to!string(f2PK.emit(B, PKptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(B, PKptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for F2 X chr male (phase-known)
  assert( f2PK.nrec(PKptg[0], PKptg[0], is_X_chr, is_female, cross_direction) == 0.0 );
  assert( f2PK.nrec(PKptg[0], PKptg[1], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(PKptg[1], PKptg[0], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(PKptg[1], PKptg[1], is_X_chr, is_female, cross_direction) == 0.0 );
}

unittest {
  writeln("Unit test " ~ __FILE__ ~ " F2 X chr female forw");

  // unit test allTrueGeno for F2 X chr female forw
  auto f2 = form_cross("F2");
  auto f2PK = form_cross_phaseknown(f2);
  auto g = f2.all_true_geno_X;
  auto gPK = f2PK.all_true_geno_X;

  assert(g.length == 3);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
  assert(g[2] == new TrueGenotype(1,1));

  // unit test allTrueGeno for F2 X chr female forw (phase-known)
  assert(gPK.length == 4);
  assert(gPK[0] == new TrueGenotype(0,0));
  assert(gPK[1] == new TrueGenotype(1,0));
  assert(gPK[2] == new TrueGenotype(0,1));
  assert(gPK[3] == new TrueGenotype(1,1));

  // X chr male
  auto is_X_chr = true;
  auto is_female = true;
  auto cross_direction = [0];

  // index to possible genotype
  auto ptgi = f2.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto ptg = new TrueGenotype[](ptgi.length);
  foreach(i, ptgival; ptgi)
    ptg[i] = g[ptgival];
  assert(ptgi.length == 2);
  assert(ptgi[0] == 0);
  assert(ptgi[1] == 1);
  assert(ptg.length == 2);
  assert(ptg[0] == new TrueGenotype(0,0));
  assert(ptg[1] == new TrueGenotype(1,0));

  auto PKptgi = f2PK.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto PKptg = new TrueGenotype[](PKptgi.length);
  foreach(i, PKptgival; PKptgi)
    PKptg[i] = gPK[PKptgival];
  assert(PKptgi.length == 2);
  assert(PKptgi[0] == 0);
  assert(PKptgi[1] == 1);
  assert(PKptg.length == 2);
  assert(PKptg[0] == new TrueGenotype(0,0));
  assert(PKptg[1] == new TrueGenotype(1,0));

  // unit test init for F2 X chr female forw
  foreach(i, gv; ptg)
    assert(to!string(f2.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  // unit test init for F2 X chr female forw (phase-known)
  foreach(gv; PKptg)
    assert(to!string(f2PK.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  //  unit test step for F2 X chr female forw
  double rf = 0.01;

  assert( to!string( f2.step(ptg[0], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));
  assert( to!string( f2.step(ptg[0], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( f2.step(ptg[1], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( f2.step(ptg[1], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));

  // unit test step for F2 X chr female forw (phase-known)
  assert( to!string( f2PK.step(PKptg[0], PKptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));
  assert( to!string( f2PK.step(PKptg[0], PKptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( f2PK.step(PKptg[1], PKptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( f2PK.step(PKptg[1], PKptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));

  // unit test emit for F2 X chr female forw
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeSymbolMapper("H");
  H ~= new TrueGenotype(1,0);
  H ~= new TrueGenotype(0,1);

  double error_prob = 0.01;

  assert(f2.emit(NA, ptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2.emit(NA, ptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2.emit(A, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(A, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(H, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(H, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));

  assert(f2PK.emit(NA, PKptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2PK.emit(NA, PKptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2PK.emit(A, PKptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(A, PKptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  assert(to!string(f2PK.emit(H, PKptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(H, PKptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for F2 X chr female forw (phase-known)
  assert( f2PK.nrec(PKptg[0], PKptg[0], is_X_chr, is_female, cross_direction) == 0.0 );
  assert( f2PK.nrec(PKptg[0], PKptg[1], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(PKptg[1], PKptg[0], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(PKptg[1], PKptg[1], is_X_chr, is_female, cross_direction) == 0.0 );
}

unittest {
  writeln("Unit test " ~ __FILE__ ~ " F2 X chr female backw");

  // unit test allTrueGeno for F2 X chr female backw
  auto f2 = form_cross("F2");
  auto f2PK = form_cross_phaseknown(f2);
  auto g = f2.all_true_geno_X;
  auto gPK = f2PK.all_true_geno_X;

  assert(g.length == 3);
  assert(g[0] == new TrueGenotype(0,0));
  assert(g[1] == new TrueGenotype(1,0));
  assert(g[2] == new TrueGenotype(1,1));

  // unit test allTrueGeno for F2 X chr female backw (phase-known)
  assert(gPK.length == 4);
  assert(gPK[0] == new TrueGenotype(0,0));
  assert(gPK[1] == new TrueGenotype(1,0));
  assert(gPK[2] == new TrueGenotype(0,1));
  assert(gPK[3] == new TrueGenotype(1,1));

  // X chr male
  auto is_X_chr = true;
  auto is_female = true;
  auto cross_direction = [1];

  // index to possible genotype
  auto ptgi = f2.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto ptg = new TrueGenotype[](ptgi.length);
  foreach(i, ptgival; ptgi)
    ptg[i] = g[ptgival];
  assert(ptgi.length == 2);
  assert(ptgi[0] == 1);
  assert(ptgi[1] == 2);
  assert(ptg.length == 2);
  assert(ptg[0] == new TrueGenotype(1,0));
  assert(ptg[1] == new TrueGenotype(1,1));

  auto PKptgi = f2PK.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto PKptg = new TrueGenotype[](PKptgi.length);
  foreach(i, PKptgival; PKptgi)
    PKptg[i] = gPK[PKptgival];
  assert(PKptgi.length == 2);
  assert(PKptgi[0] == 2);
  assert(PKptgi[1] == 3);
  assert(PKptg.length == 2);
  assert(PKptg[0] == new TrueGenotype(0,1));
  assert(PKptg[1] == new TrueGenotype(1,1));

  // unit test init for F2 X chr female backw
  foreach(i, gv; ptg)
    assert(to!string(f2.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  // unit test init for F2 X chr female backw (phase-known)
  foreach(gv; PKptg)
    assert(to!string(f2PK.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  //  unit test step for F2 X chr female backw
  double rf = 0.01;

  assert( to!string( f2.step(ptg[0], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));
  assert( to!string( f2.step(ptg[0], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( f2.step(ptg[1], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( f2.step(ptg[1], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));

  // unit test step for F2 X chr female backw (phase-known)
  assert( to!string( f2PK.step(PKptg[0], PKptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));
  assert( to!string( f2PK.step(PKptg[0], PKptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));

  assert( to!string( f2PK.step(PKptg[1], PKptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(rf) ));
  assert( to!string( f2PK.step(PKptg[1], PKptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(1-rf) ));

  // unit test emit for F2 X chr female backw
  auto NA = new GenotypeSymbolMapper("-");
  auto B = new GenotypeSymbolMapper("B");
  B ~= new TrueGenotype(1,1);
  auto H = new GenotypeSymbolMapper("H");
  H ~= new TrueGenotype(1,0);
  H ~= new TrueGenotype(0,1);

  double error_prob = 0.01;

  assert(f2.emit(NA, ptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2.emit(NA, ptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2.emit(H, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(H, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(B, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(B, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));

  assert(f2PK.emit(NA, PKptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(f2PK.emit(NA, PKptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);

  assert(to!string(f2PK.emit(H, PKptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(H, PKptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  assert(to!string(f2PK.emit(B, PKptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(B, PKptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for F2 X chr female backw (phase-known)
  assert( f2PK.nrec(PKptg[0], PKptg[0], is_X_chr, is_female, cross_direction) == 0.0 );
  assert( f2PK.nrec(PKptg[0], PKptg[1], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(PKptg[1], PKptg[0], is_X_chr, is_female, cross_direction) == 1.0 );
  assert( f2PK.nrec(PKptg[1], PKptg[1], is_X_chr, is_female, cross_direction) == 0.0 );
}

unittest {
  writeln("Unit test " ~ __FILE__ ~ " RISIB autosome");

  // unit test all_true_geno for RISIB, autosome
  auto ril = form_cross("RISIB");
  auto atg = ril.all_true_geno_A;
  assert(atg.length == 2);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,1));

  // blank stuff for autosome
  auto is_X_chr = false;
  auto is_female = false;
  auto cross_direction = new int[](1);

  // unit test init for RISIB autosome
  foreach(gv; atg)
    assert(to!string(ril.init(gv, is_X_chr, is_female, cross_direction)) == to!string(log(0.5)));

  // unit test step for RISIB autosome
  double rf = 0.01, R = 4.0*rf/(1.0+6.0*rf);

  assert( to!string( ril.step(atg[0], atg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-R)) ));
  assert( to!string( ril.step(atg[0], atg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(R) ));

  assert( to!string( ril.step(atg[1], atg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(R) ));
  assert( to!string( ril.step(atg[1], atg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log((1-R)) ));

  // unit test emit for RISIB autosome
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto B = new GenotypeSymbolMapper("B");
  B ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(ril.emit(NA, atg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(ril.emit(NA, atg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(to!string(ril.emit(A, atg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(ril.emit(A, atg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(ril.emit(B, atg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(ril.emit(B, atg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for RISIB autosome
  assert(ril.nrec(atg[0], atg[0], is_X_chr, is_female, cross_direction) == 0.0);
  assert(ril.nrec(atg[0], atg[1], is_X_chr, is_female, cross_direction) == 1.0);

  assert(ril.nrec(atg[1], atg[0], is_X_chr, is_female, cross_direction) == 1.0);
  assert(ril.nrec(atg[1], atg[1], is_X_chr, is_female, cross_direction) == 0.0);
}

unittest {
  writeln("Unit test " ~ __FILE__ ~ " RISIB X chr forward");

  // unit test all_true_geno for RISIB, X chr, forward
  auto ril = form_cross("RISIB");
  auto atg = ril.all_true_geno_X;
  assert(atg.length == 2);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,1));

  // X chr forward
  auto is_X_chr = true;
  auto is_female = false;
  auto cross_direction = [0];

  // index to possible genotype
  auto ptgi = ril.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto ptg = new TrueGenotype[](ptgi.length);
  foreach(i, ptgival; ptgi)
    ptg[i] = atg[ptgival];
  assert(ptgi.length == 2);
  assert(ptgi[0] == 0);
  assert(ptgi[1] == 1);
  assert(ptg.length == 2);
  assert(ptg[0] == new TrueGenotype(0,0));
  assert(ptg[1] == new TrueGenotype(1,1));

  // unit test init for RISIB X chr, forward
  assert(to!string(ril.init(ptg[0], is_X_chr, is_female, cross_direction)) == to!string(log(2.0/3.0)));
  assert(to!string(ril.init(ptg[1], is_X_chr, is_female, cross_direction)) == to!string(log(1.0/3.0)));

  // unit test step for RISIB X chr, forward
  double rf = 0.01;
  double P00 = (1.0+2.0*rf)/(1.0+4.0*rf), P01 = 2.0*rf/(1.0+4.0*rf),
    P10 = 4.0*rf/(1.0+4.0*rf), P11 = 1.0/(1.0+4.0*rf);

  assert( to!string( ril.step(ptg[0], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P00) ));
  assert( to!string( ril.step(ptg[0], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P01) ));

  assert( to!string( ril.step(ptg[1], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P10) ));
  assert( to!string( ril.step(ptg[1], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P11) ));

  // unit test emit for RISIB X chr, forward
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto B = new GenotypeSymbolMapper("B");
  B ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(ril.emit(NA, ptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(ril.emit(NA, ptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(to!string(ril.emit(A, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(ril.emit(A, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(ril.emit(B, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(ril.emit(B, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for RISIB X chr, forward
  assert(ril.nrec(ptg[0], ptg[0], is_X_chr, is_female, cross_direction) == 0.0);
  assert(ril.nrec(ptg[0], ptg[1], is_X_chr, is_female, cross_direction) == 1.0);

  assert(ril.nrec(ptg[1], ptg[0], is_X_chr, is_female, cross_direction) == 1.0);
  assert(ril.nrec(ptg[1], ptg[1], is_X_chr, is_female, cross_direction) == 0.0);
}

unittest {
  writeln("Unit test " ~ __FILE__ ~ " RISIB X chr reverse");

  // unit test all_true_geno for RISIB, X chr, reverse
  auto ril = form_cross("RISIB");
  auto atg = ril.all_true_geno_X;
  assert(atg.length == 2);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,1));

  // X chr reverse
  auto is_X_chr = true;
  auto is_female = false;
  auto cross_direction = [1];

  // index to possible genotype
  auto ptgi = ril.possible_true_geno_index(is_X_chr, is_female, cross_direction);
  auto ptg = new TrueGenotype[](ptgi.length);
  foreach(i, ptgival; ptgi)
    ptg[i] = atg[ptgival];
  assert(ptgi.length == 2);
  assert(ptgi[0] == 0);
  assert(ptgi[1] == 1);
  assert(ptg.length == 2);
  assert(ptg[0] == new TrueGenotype(0,0));
  assert(ptg[1] == new TrueGenotype(1,1));

  // unit test init for RISIB X chr, reverse
  assert(to!string(ril.init(ptg[0], is_X_chr, is_female, cross_direction)) == to!string(log(1.0/3.0)));
  assert(to!string(ril.init(ptg[1], is_X_chr, is_female, cross_direction)) == to!string(log(2.0/3.0)));

  // unit test step for RISIB X chr, reverse
  double rf = 0.01;
  double P11 = (1.0+2.0*rf)/(1.0+4.0*rf), P10 = 2.0*rf/(1.0+4.0*rf),
    P01 = 4.0*rf/(1.0+4.0*rf), P00 = 1.0/(1.0+4.0*rf);

  assert( to!string( ril.step(ptg[0], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P00) ));
  assert( to!string( ril.step(ptg[0], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P01) ));

  assert( to!string( ril.step(ptg[1], ptg[0], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P10) ));
  assert( to!string( ril.step(ptg[1], ptg[1], rf, is_X_chr, is_female, cross_direction) ) ==
          to!string( log(P11) ));

  // unit test emit for RISIB X chr, reverse
  auto NA = new GenotypeSymbolMapper("-");
  auto A = new GenotypeSymbolMapper("A");
  A ~= new TrueGenotype(0,0);
  auto B = new GenotypeSymbolMapper("B");
  B ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(ril.emit(NA, ptg[0], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(ril.emit(NA, ptg[1], error_prob, is_X_chr, is_female, cross_direction) == 0);
  assert(to!string(ril.emit(A, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(ril.emit(A, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));
  assert(to!string(ril.emit(B, ptg[1], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(1.0-error_prob)));
  assert(to!string(ril.emit(B, ptg[0], error_prob, is_X_chr, is_female, cross_direction)) == to!string(log(error_prob)));

  // unit test nrec for RISIB X chr, reverse
  assert(ril.nrec(ptg[0], ptg[0], is_X_chr, is_female, cross_direction) == 0.0);
  assert(ril.nrec(ptg[0], ptg[1], is_X_chr, is_female, cross_direction) == 1.0);

  assert(ril.nrec(ptg[1], ptg[0], is_X_chr, is_female, cross_direction) == 1.0);
  assert(ril.nrec(ptg[1], ptg[1], is_X_chr, is_female, cross_direction) == 0.0);
}
