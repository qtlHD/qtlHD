/**
 * HMM module
 */

module qtl.core.hmm;

// things I think I really need
import qtl.core.primitives;
import qtl.core.genotype;
import std.stdio;
import std.math;
import qtl.core.hmm_f2;

// things for the unit tests 
import qtl.plugins.input.read_csv;
import std.path;


// calculate QTL genotype probabilities
double[int][int][T] calcGenoprob(T)(T[int][int] genotypes, double[] rec_frac, double error_prob)
in {
  assert(genotypes.length == rec_frac.length-1);
 }
body {
  double[int][T] alpha, beta;
  double[int][int][T] genoprobs;

  auto n_markers = geno[0].length;

  foreach(i; 0..geno.length) {

    // initialize alpha and beta
    //    foreach(true_gen; cross.possible_true_genotypes) {
    //      alpha[true_gen][0] = cross.init(true_gen) + cross.emit(geno[i][0], true_gen, error_prob);
    //      beta[true_gen][n_marekrs-1] = 0.0;
    //    }
  }
  return genoprobs;
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","test","data","input","listeria.csv");
  writeln("  - read CSV " ~ fn);
  auto data = new ReadSimpleCSV!F2(fn);
  //  auto cross = new F2Cross(data.genotypes);

  //  auto genoprob = calcGenoprob!F2(cross, 0, 0.01);
}



