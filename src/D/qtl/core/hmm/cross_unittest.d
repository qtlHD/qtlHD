/**
 * Unittests for cross.d
 **/

module qtl.core.hmm.cross_unittest;

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
  auto atg = bc.all_true_geno(false);
  assert(atg.length == 2);
  assert(atg[0] == new TrueGenotype(0,0));
  assert(atg[1] == new TrueGenotype(1,0));

  auto is_X_chr = false;
  auto is_female = false;
  auto cross_direction = new int[][](1,1);

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
  auto NA = new GenotypeCombinator("-");
  auto A = new GenotypeCombinator("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeCombinator("H");
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
  writeln("Unit test " ~ __FILE__ ~ " F2");

  // unit test allTrueGeno for F2
  auto f2 = form_cross("F2");
  auto f2PK = form_cross_phaseknown(f2);
  auto g = f2.all_true_geno;
  auto gPK = f2PK.all_true_geno;

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

  // unit test init for F2
  double val[] = [log(0.25), log(0.5), log(0.25)];

  foreach(i, gv; g)
    assert(to!string(f2.init(gv)) == to!string(val[i]));

  // unit test init for F2 (phase-known)
  foreach(gv; gPK)
    assert(to!string(f2PK.init(gv)) == to!string(val[0]));

  //  unit test step for F2
  double rf = 0.01;

  assert( to!string( f2.step(g[0], g[0], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2.step(g[0], g[1], rf) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( f2.step(g[0], g[2], rf) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( f2.step(g[1], g[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2.step(g[1], g[1], rf) ) ==
          to!string( log(rf^^2 + (1-rf)^^2) ));
  assert( to!string( f2.step(g[1], g[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));

  assert( to!string( f2.step(g[2], g[0], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2.step(g[2], g[1], rf) ) ==
          to!string( log(2.0*rf*(1-rf)) ));
  assert( to!string( f2.step(g[2], g[2], rf) ) ==
          to!string( log((1-rf)^^2) ));


  // unit test step for F2 (phase-known)
  assert( to!string( f2PK.step(gPK[0], gPK[0], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[0], gPK[1], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[0], gPK[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[0], gPK[3], rf) ) ==
          to!string( log(rf^^2) ));

  assert( to!string( f2PK.step(gPK[1], gPK[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[1], gPK[1], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[1], gPK[2], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2PK.step(gPK[1], gPK[3], rf) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( f2PK.step(gPK[2], gPK[0], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[2], gPK[2], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[2], gPK[1], rf) ) ==
          to!string( log(rf^^2) ));
  assert( to!string( f2PK.step(gPK[2], gPK[3], rf) ) ==
          to!string( log(rf*(1.0-rf)) ));

  assert( to!string( f2PK.step(gPK[3], gPK[3], rf) ) ==
          to!string( log((1-rf)^^2) ));
  assert( to!string( f2PK.step(gPK[3], gPK[1], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[3], gPK[2], rf) ) ==
          to!string( log(rf*(1-rf)) ));
  assert( to!string( f2PK.step(gPK[3], gPK[0], rf) ) ==
          to!string( log(rf^^2) ));

  // unit test emit for F2
  auto NA = new GenotypeCombinator("-");
  auto A = new GenotypeCombinator("A");
  A ~= new TrueGenotype(0,0);
  auto H = new GenotypeCombinator("H");
  H ~= new TrueGenotype(1,0);
  H ~= new TrueGenotype(0,1);
  auto B = new GenotypeCombinator("B");
  B ~= new TrueGenotype(1,1);
  auto AorH = new GenotypeCombinator("AorH");
  AorH ~= new TrueGenotype(0,0);
  AorH ~= new TrueGenotype(1,0);
  AorH ~= new TrueGenotype(0,1);
  auto HorB = new GenotypeCombinator("HorB");
  HorB ~= new TrueGenotype(1,0);
  HorB ~= new TrueGenotype(0,1);
  HorB ~= new TrueGenotype(1,1);

  double error_prob = 0.01;

  assert(f2.emit(NA, g[0], error_prob) == 0);
  assert(f2.emit(NA, g[1], error_prob) == 0);
  assert(f2.emit(NA, g[2], error_prob) == 0);

  assert(to!string(f2.emit(A, g[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(A, g[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(A, g[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(H, g[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(H, g[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(H, g[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(B, g[2], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2.emit(B, g[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2.emit(B, g[1], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2.emit(AorH, g[0], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, g[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(AorH, g[2], error_prob)) == to!string(log(error_prob)));

  assert(to!string(f2.emit(HorB, g[0], error_prob)) == to!string(log(error_prob)));
  assert(to!string(f2.emit(HorB, g[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2.emit(HorB, g[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));

  assert(f2PK.emit(NA, gPK[0], error_prob) == 0);
  assert(f2PK.emit(NA, gPK[1], error_prob) == 0);
  assert(f2PK.emit(NA, gPK[2], error_prob) == 0);
  assert(f2PK.emit(NA, gPK[3], error_prob) == 0);

  assert(to!string(f2PK.emit(A, gPK[0], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(A, gPK[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(A, gPK[2], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(A, gPK[3], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2PK.emit(H, gPK[1], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(H, gPK[2], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(H, gPK[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(H, gPK[3], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2PK.emit(B, gPK[3], error_prob)) == to!string(log(1.0-error_prob)));
  assert(to!string(f2PK.emit(B, gPK[0], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(B, gPK[1], error_prob)) == to!string(log(error_prob/2.0)));
  assert(to!string(f2PK.emit(B, gPK[2], error_prob)) == to!string(log(error_prob/2.0)));

  assert(to!string(f2PK.emit(AorH, gPK[0], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(AorH, gPK[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(AorH, gPK[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(AorH, gPK[3], error_prob)) == to!string(log(error_prob)));

  assert(to!string(f2PK.emit(HorB, gPK[0], error_prob)) == to!string(log(error_prob)));
  assert(to!string(f2PK.emit(HorB, gPK[1], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(HorB, gPK[2], error_prob)) == to!string(log(1.0-error_prob/2.0)));
  assert(to!string(f2PK.emit(HorB, gPK[3], error_prob)) == to!string(log(1.0-error_prob/2.0)));

  // unit test nrec for F2 (phase-known)
  assert( f2PK.nrec(gPK[0], gPK[0]) == 0.0 );
  assert( f2PK.nrec(gPK[0], gPK[1]) == 0.5 );
  assert( f2PK.nrec(gPK[0], gPK[2]) == 0.5 );
  assert( f2PK.nrec(gPK[0], gPK[3]) == 1.0 );

  assert( f2PK.nrec(gPK[1], gPK[0]) == 0.5 );
  assert( f2PK.nrec(gPK[1], gPK[1]) == 0.0 );
  assert( f2PK.nrec(gPK[1], gPK[2]) == 1.0 );
  assert( f2PK.nrec(gPK[1], gPK[3]) == 0.5 );

  assert( f2PK.nrec(gPK[2], gPK[0]) == 0.5 );
  assert( f2PK.nrec(gPK[2], gPK[2]) == 0.0 );
  assert( f2PK.nrec(gPK[2], gPK[1]) == 1.0 );
  assert( f2PK.nrec(gPK[2], gPK[3]) == 0.5 );

  assert( f2PK.nrec(gPK[3], gPK[3]) == 0.0 );
  assert( f2PK.nrec(gPK[3], gPK[1]) == 0.5 );
  assert( f2PK.nrec(gPK[3], gPK[2]) == 0.5 );
  assert( f2PK.nrec(gPK[3], gPK[0]) == 1.0 );
}
