/**
 * calc_genoprob
 */

module qtl.core.hmm_calcgenoprob;

import std.string;
import std.conv;

// calculate QTL genotype probabilities
mixin template calcGenoprobCode(T)
{
  double[T][int][int] calc_geno_prob(Genotype!T[][] genotypes, double[] rec_frac, double error_prob)
  {
    if(genotypes[0].length != rec_frac.length+1) {
      // throw new Exception("no. markers in genotypes " ~ to!string(genotypes[0].length) ~ "doesn't match rec_frac length" ~ to!string(rec_frac.length+1));
      writeln(genotypes[0].length);
      writeln(rec_frac.length+1);
      throw new Exception("no. markers in genotypes doesn't match rec_frac length");
    }
    if(error_prob < 0.0 || error_prob > 1.0)
      throw new Exception("error_prob out of range");
    foreach(rf; rec_frac) {
      if(rf < 0 || rf > 0.5)
        throw new Exception("rec_frac must be >= 0 and <= 0.5");
    }

    int n_individuals = genotypes.length;
    int n_markers = genotypes[0].length;
    auto all_true_geno = allTrueGeno(genotypes[0][0].value);

    double[int][T] alpha, beta;
    double[T][int][int] genoprobs;

    foreach(ind; 0..n_individuals) {
      alpha = forwardEquations(genotypes[ind], all_true_geno, rec_frac, error_prob);
        
      beta = backwardEquations(genotypes[ind], all_true_geno, rec_frac, error_prob);

      // calculate genotype probabilities
      double sum_at_pos;
      foreach(pos; 0..n_markers) {
        sum_at_pos = genoprobs[ind][pos][all_true_geno[0]] = alpha[all_true_geno[0]][pos] + beta[all_true_geno[0]][pos];
        foreach(true_geno; all_true_geno[1..$]) {
          genoprobs[ind][pos][true_geno] = alpha[true_geno][pos] + beta[true_geno][pos];
          sum_at_pos = addlog(sum_at_pos, genoprobs[ind][pos][true_geno]);
        }
        foreach(true_geno; all_true_geno) {
          genoprobs[ind][pos][true_geno] = exp(genoprobs[ind][pos][true_geno] - sum_at_pos);
        }
      }
    }

    return genoprobs;
  }
}

