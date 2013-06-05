/**
 * Cross classes for HMM-related functions
 **/

module qtl.core.hmm.cross;

import std.math, std.stdio, std.path;
import std.exception;
import std.conv;

import qtl.core.primitives;
import qtl.core.genotype;

// class to contain HMM-related functions
class Cross {
  string cross_type, phase_known_class_name;

  TrueGenotype[] all_true_geno_A, all_true_geno_X;

  abstract TrueGenotype[] all_true_geno(bool is_X_chr);

  abstract double init(T)(TrueGenotype truegen, bool is_X_chr, bool is_female, T cross_direction);
  abstract double step(T)(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac,
                          bool is_X_chr, bool is_female, T cross_direction);
  abstract double emit(T)(GenotypeCombinator obsgen, TrueGenotype truegen, double error_prob,
                          bool is_X_chr, bool is_female, T cross_direction);
  abstract double nrec(T)(TrueGenotype truegen_left, TrueGenotype truegen_right,
                          bool is_X_chr, bool is_female, T cross_direction);

  TrueGenotype[] all_true_geno(bool is_X_chr)
  {
    if(is_X_chr) return all_true_geno_X;
    else return all_true_geno_A;
  }

  abstract size_t[] possible_true_geno_index(T)(bool is_X_chr, bool is_female, T cross_direction);
}

// the switch
Cross form_cross(string which_cross) {

  Cross cross = cast(Cross)Object.factory("qtl.core.hmm.cross." ~ which_cross);

  if(!cross) {
    throw new Exception("cross type " ~ which_cross ~ " not supported.");
  }

  return cross;
}

// the second switch: phase-known versions
Cross form_cross_phaseknown(Cross cross) {

  string which_cross = cross.phase_known_class_name;
  if(!which_cross) {
    throw new Exception("input cross doesn't contain phase_known_class_name\n" ~
                        "    to indicate corresponding phase-known cross class.");
  }

  Cross pkcross = cast(Cross)Object.factory("qtl.core.hmm.cross." ~ which_cross);

  if(!pkcross) {
    throw new Exception("cross type " ~ cross.cross_type ~ " not supported.");
  }

  return(pkcross);
}

// BC (backcross)
class BC : Cross {
  this() {
    cross_type = "BC";

    // vector of true genotypes
    auto g1 = new TrueGenotype(0,0);
    auto g2 = new TrueGenotype(1,0);
    auto g3 = new TrueGenotype(1,1);
    all_true_geno_A = [g1, g2];
    all_true_geno_X = [g1, g2, g3];

    phase_known_class_name = "BC";
  }

  // indexes to possible true genotypes, by chromosome type and sex
  size_t[] possible_true_geno_index(bool is_X_chr, bool is_female, bool ignored_argument)
  {
    if(is_X_chr) {
      if(is_female) return([0,1]);
      else return([0,2]);
    }
    else {
      return([0, 1]);
    }
  }



  // ln Pr(true genotype)
  double init(TrueGenotype truegen, bool is_X_chr, bool is_female, bool ignored_argument)
  {
    if(is_X_chr) {
      if(is_female) return(initA(truegen));
      else return(initXmale(truegen));
    }
    else {
      return initA(truegen);
    }
  }

  double initA(TrueGenotype truegen)
  {
    if(truegen != all_true_geno_A[0] &&
       truegen != all_true_geno_A[1])
      throw new Exception("truegen must be 0,0 or 1,0");

    return(-LN2);
  }

  double initXmale(TrueGenotype truegen)
  {
    if(truegen != all_true_geno_X[0] &&
       truegen != all_true_geno_X[2])
      throw new Exception("Male on X: truegen must be 0,0 or 1,1");

    return(-LN2);
  }



  // ln Pr(genotype at right marker | genotype at left marker)
  double step(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac,
              bool is_X_chr, bool is_female, bool ignored_argument)
  {
    if(is_X_chr) {
      if(is_female) return(stepA(truegen_left, truegen_right, rec_frac));
      else return(stepXmale(truegen_left, truegen_right, rec_frac));
    }
    else {
      return(stepA(truegen_left, truegen_right, rec_frac));
    }
  }

  double stepA(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_A atg;

    if(truegen_left == atg[0]) {
      if(truegen_right == atg[0]) // A -> A
        return(log(1.0-rec_frac));
      if(truegen_right == atg[1]) // A -> H
        return(log(rec_frac));
    }

    if(truegen_left == atg[1]) {
      if(truegen_right == atg[1]) // H -> H
        return(log(1.0-rec_frac));
      if(truegen_right == atg[0]) // H -> A
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepXmale(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_X atg;

    if(truegen_left == atg[0]) {
      if(truegen_right == atg[0]) // A -> A
        return(log(1.0-rec_frac));
      if(truegen_right == atg[2]) // A -> B
        return(log(rec_frac));
    }

    if(truegen_left == atg[2]) {
      if(truegen_right == atg[2]) // B -> B
        return(log(1.0-rec_frac));
      if(truegen_right == atg[0]) // B -> A
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }



  // ln Pr(observed genotype | true genotype)
  double emit(GenotypeCombinator obsgen, TrueGenotype truegen, double error_prob,
              bool is_X_chr, bool is_female, bool ignored_argument) // don't actually use X chr or sex here
  {
    if(obsgen.list.length==0) // missing value
      return(0.0); // log(1.0)

    if(obsgen.match(truegen)) // compatible with truegen?
      return(log(1.0-error_prob));
    else
      return(log(error_prob));
  }



  // No. recombination events
  double nrec(TrueGenotype truegen_left, TrueGenotype truegen_right,
              bool is_X_chr, bool is_female, bool ignored_argument)
  {
    if(is_X_chr) {
      if(is_female) return(nrecA(truegen_left, truegen_right));
      else return(nrecXmale(truegen_left, truegen_right));
    }
    else {
      return(nrecA(truegen_left, truegen_right, rec_frac));
    }
  }

  double nrecA(TrueGenotype truegen_left, TrueGenotype truegen_right)
  {
    alias all_true_geno_A atg;

    if(truegen_left == atg[0]) {
      if(truegen_right == atg[0]) // A -> A
        return(0.0);
      if(truegen_right == atg[1]) // A -> H
        return(1.0);
    }

    if(truegen_left == atg[1]) {
      if(truegen_right == atg[1]) // H -> H
        return(0.0);
      if(truegen_right == atg[0]) // H -> A
        return(1.0);
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double nrecXmale(TrueGenotype truegen_left, TrueGenotype truegen_right)
  {
    alias all_true_geno_X atg;

    if(truegen_left == atg[0]) {
      if(truegen_right == atg[0]) // A -> A
        return(0.0);
      if(truegen_right == atg[2]) // A -> B
        return(1.0);
    }

    if(truegen_left == atg[2]) {
      if(truegen_right == atg[2]) // B -> B
        return(0.0);
      if(truegen_right == atg[0]) // B -> A
        return(1.0);
    }

    throw new Exception("inputs not among the possible true genotypes");
  }
}



// F2 (intercross) [ used in calc_geno_prob]
class F2 : Cross {

  this() {
    cross_type = "F2";

    auto g1 = new TrueGenotype(0,0);
    auto g2 = new TrueGenotype(1,0);
    auto g3 = new TrueGenotype(1,1);
    all_true_geno_X = all_true_geno_A = [g1, g2, g3];

    phase_known_class_name = "F2_phaseknown";
  }

  // indexes to possible true genotypes, by chromosome type and sex
  size_t[] possible_true_geno_index(bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    if(is_X_chr) {
      if(is_female) {
        if(is_forward_cross) return([0,1]);
        else return([1,2]);
      }
      else return([0,2]);
    }
    else {
      return([0,1,2]);
    }
  }

  // ln Pr(true genotype)
  double init(TrueGenotype truegen, bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    if(is_X_chr) {
      if(is_female)
        return(initXfemale(truegen, is_forward_cross));
      else
        return(initXmale(truegen));
    }
    else return(initA(truegen));
  }

  double initA(TrueGenotype truegen)
  {
    alias all_true_geno_A atg;

    if(truegen==atg[0] || truegen==atg[2]) // AA or BB
      return(-2.0*LN2);

    if(truegen==atg[1]) // AB
      return(-LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  double initXmale(TrueGenotype truegen)
  {
    alias all_true_geno_X atg;

    if(truegen==atg[0] || truegen==atg[2]) // A- or B-
      return(-LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  double initXfemale(TrueGenotype truegen, bool is_forward_cross)
  {
    alias all_true_geno_X atg;

    if((is_forward_cross && (truegen==atg[0] || truegen==atg[1])) ||
       (!is_forward_cross && (truegen==atg[1] || truegen==atg[2])))
      return(-LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  // ln Pr(genotype at right marker | genotype at left marker)
  double step(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac,
              bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    if(is_X_chr) {
      if(is_female) {
        if(is_forward_cross) return(stepXfemaleforw(truegen_left, truegen_right, rec_frac));
        else return(stepXfemalerev(truegen_left, truegen_right, rec_frac));
      }
      else return(stepXmale(truegen_left, truegen_right, rec_frac));
    }
    else return(stepA(truegen_left, truegen_right, rec_frac));
  }

  double stepA(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_A atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // A -> A
        return(2.0*log(1.0-rec_frac));
      if(truegen_right==atg[1]) // A -> H
        return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atg[2]) // A -> B
        return(2.0*log(rec_frac));
    }
    else if(truegen_left == atg[1]) {
      if(truegen_right==atg[0] || // H -> A
         truegen_right==atg[2]) // H -> B
        return(log(rec_frac) + log(1.0-rec_frac));
      if(truegen_right==atg[1]) // H -> H
        return(log((1.0-rec_frac)^^2 + rec_frac^^2));
    }
    else if(truegen_left == atg[2]) {
      if(truegen_right==atg[0]) // B -> A
        return(2.0*log(rec_frac));
      if(truegen_right==atg[1]) // B -> H
        return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atg[2]) // B -> B
        return(2.0*log(1.0-rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepA(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_A atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // A -> A
        return(2.0*log(1.0-rec_frac));
      if(truegen_right==atg[1]) // A -> H
        return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atg[2]) // A -> B
        return(2.0*log(rec_frac));
    }
    else if(truegen_left == atg[1]) {
      if(truegen_right==atg[0] || // H -> A
         truegen_right==atg[2]) // H -> B
        return(log(rec_frac) + log(1.0-rec_frac));
      if(truegen_right==atg[1]) // H -> H
        return(log((1.0-rec_frac)^^2 + rec_frac^^2));
    }
    else if(truegen_left == atg[2]) {
      if(truegen_right==atg[0]) // B -> A
        return(2.0*log(rec_frac));
      if(truegen_right==atg[1]) // B -> H
        return(LN2 + log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atg[2]) // B -> B
        return(2.0*log(1.0-rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepXmale(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // AA -> AA
        return(log(1.0-rec_frac));
      if(truegen_right==atg[2]) // AA -> BB
        return(log(rec_frac));
    }
    else if(truegen_left == atg[2]) {
      if(truegen_right==atg[2]) // BB -> BB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[0]) // BB -> AA
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepXfemaleforw(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // AA -> AA
        return(log(1.0-rec_frac));
      if(truegen_right==atg[1]) // AA -> AB
        return(log(rec_frac));
    }
    else if(truegen_left == atg[1]) {
      if(truegen_right==atg[1]) // AB -> AB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[0]) // AB -> AA
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepXfemalerev(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[1]) {
      if(truegen_right==atg[1]) // AB -> AB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[2]) // AB -> BB
        return(log(rec_frac));
    }
    else if(truegen_left == atg[2]) {
      if(truegen_right==atg[2]) // BB -> BB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[1]) // BB -> AB
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  // ln Pr(observed genotype | true genotype)
  double emit(GenotypeCombinator obsgen, TrueGenotype truegen, double error_prob,
              bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    auto n_obsgen = obsgen.list.length;

    if(n_obsgen==0) // missing value
      return(0.0); // log(1.0)

    if(obsgen.match(truegen)) { // compatible with truegen
      if(n_obsgen>2) // AorH and HorB cases
        return(log(1.0-error_prob/2.0));
      else // A, H, B cases
        return(log(1.0-error_prob));
    }
    else {
      if(n_obsgen>2) // AorH and HorB cases
        return(log(error_prob));
      else // A, H, B cases
        return(log(error_prob/2.0));
    }
  }

  // proportion of recombination events
  double nrec(TrueGenotype truegen_left, TrueGenotype truegen_right,
              bool is_X_chr, bool is_female, bool forw_direction)
  {
    throw new Exception("nrec undefined for cross type F2");
  }
}

// F2_phaseknown (intercross with phase-known genotypes) [used in est_map]
class F2_phaseknown : F2 {
  this() {
    cross_type = "F2_phaseknown";

    auto g1 = new TrueGenotype(0,0);
    auto g2 = new TrueGenotype(1,0);
    auto g3 = new TrueGenotype(0,1);
    auto g4 = new TrueGenotype(1,1);
    all_true_geno_X = all_true_geno_A = [g1, g2, g3, g4];
  }

  // indexes to possible true genotypes, by chromosome type and sex
  size_t[] possible_true_geno_index(bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    if(is_X_chr) {
      if(is_female) {
        if(is_forward_cross) return([0,1]);
        else return([2,3]);
      }
      else return([0,3]);
    }
    else {
      return([0,1,2,3]);
    }
  }

  // ln Pr(true genotype)
  double init(TrueGenotype truegen, bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    if(is_X_chr) {
      if(is_female)
        return(initXfemale(truegen, is_forward_cross));
      else
        return(initXmale(truegen));
    }
    else return(initA(truegen));
  }

  double initA(TrueGenotype truegen)
  {
    alias all_true_geno_A atg;

    if(truegen==atg[0] || truegen==atg[1] ||
       truegen==atg[2] || truegen==atg[3])
      return(-2.0*LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  double initXmale(TrueGenotype truegen)
  {
    alias all_true_geno_X atg;

    if(truegen==atg[0] || truegen==atg[3]) // A- or B-
      return(-LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  double initXfemale(TrueGenotype truegen, bool is_forward_cross)
  {
    alias all_true_geno_X atg;

    if((is_forward_cross && (truegen==atg[0] || truegen==atg[1])) ||
       (!is_forward_cross && (truegen==atg[2] || truegen==atg[3])))
      return(-LN2);

    throw new Exception("truegen not among possible true genotypes");
  }

  // ln Pr(genotype at right marker | genotype at left marker)
  double step(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac,
              bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    if(is_X_chr) {
      if(is_female) {
        if(is_forward_cross) return(stepXfemaleforw(truegen_left, truegen_right, rec_frac));
        else return(stepXfemalerev(truegen_left, truegen_right, rec_frac));
      }
      else return(stepXmale(truegen_left, truegen_right, rec_frac));
    }
    else return(stepA(truegen_left, truegen_right, rec_frac));
  }

  double stepA(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno atgpk;

    if(truegen_left==atgpk[0]) {
      if(truegen_right==atgpk[0]) // AA -> AA
        return(2.0*log(1.0-rec_frac));
      if(truegen_right==atgpk[1] || truegen_right==atgpk[2]) // AA -> AB; AA -> BA
        return(log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atgpk[3]) // AA -> BB
        return(2.0*log(rec_frac));
    }
    else if(truegen_left == atgpk[1]) {
      if(truegen_right==atgpk[0] || // AB -> AA
         truegen_right==atgpk[3]) // AB -> BB
        return(log(rec_frac) + log(1.0-rec_frac));
      if(truegen_right==atgpk[1]) // AB -> AB
        return(log((1.0-rec_frac)^^2));
      if(truegen_right==atgpk[2]) // AB -> BA
        return(log(rec_frac^^2));
    }
    else if(truegen_left == atgpk[2]) {
      if(truegen_right==atgpk[0] || // BA -> AA
         truegen_right==atgpk[3]) // BA -> BB
        return(log(rec_frac) + log(1.0-rec_frac));
      if(truegen_right==atgpk[2]) // BA -> BA
        return(log((1.0-rec_frac)^^2));
      if(truegen_right==atgpk[1]) // BA -> AB
        return(log(rec_frac^^2));
    }
    else if(truegen_left == atgpk[3]) {
      if(truegen_right==atgpk[0]) // BB -> AA
        return(2.0*log(rec_frac));
      if(truegen_right==atgpk[1] || // BB -> AB
         truegen_right==atgpk[2])   // BB -> BA
        return(log(1.0-rec_frac) + log(rec_frac));
      if(truegen_right==atgpk[3]) // BB -> BB
        return(2.0*log(1.0-rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepXmale(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // AA -> AA
        return(log(1.0-rec_frac));
      if(truegen_right==atg[3]) // AA -> BB
        return(log(rec_frac));
    }
    else if(truegen_left == atg[3]) {
      if(truegen_right==atg[3]) // BB -> BB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[0]) // BB -> AA
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepXfemaleforw(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // AA -> AA
        return(log(1.0-rec_frac));
      if(truegen_right==atg[1]) // AA -> AB
        return(log(rec_frac));
    }
    else if(truegen_left == atg[1]) {
      if(truegen_right==atg[1]) // AB -> AB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[0]) // AB -> AA
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double stepXfemalerev(TrueGenotype truegen_left, TrueGenotype truegen_right, double rec_frac)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[2]) {
      if(truegen_right==atg[2]) // AB -> AB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[3]) // AB -> BB
        return(log(rec_frac));
    }
    else if(truegen_left == atg[3]) {
      if(truegen_right==atg[3]) // BB -> BB
        return(log(1.0-rec_frac));
      if(truegen_right==atg[2]) // BB -> AB
        return(log(rec_frac));
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  // proportion of recombination events
  double nrec(TrueGenotype truegen_left, TrueGenotype truegen_right,
              bool is_X_chr, bool is_female, bool is_forward_cross)
  {
    if(is_X_chr) {
      if(is_female) {
        if(is_forward_cross) return(nrecXfemaleforw(truegen_left, truegen_right));
        else return(nrecXfemalerev(truegen_left, truegen_right));
      }
      else return(nrecXmale(truegen_left, truegen_right));
    }
    else return(nrecA(truegen_left, truegen_right));
  }

  double nrecA(TrueGenotype truegen_left, TrueGenotype truegen_right)
  {
    alias all_true_geno_A atgpk;

    if(truegen_left==atgpk[0]) {
      if(truegen_right==atgpk[0]) // AA -> AA
        return(0.0);
      if(truegen_right==atgpk[1] || truegen_right==atgpk[2]) // AA -> AB; AA -> BA
        return(0.5);
      if(truegen_right==atgpk[3]) // AA -> BB
        return(1.0);
    }
    else if(truegen_left == atgpk[1]) {
      if(truegen_right==atgpk[0] || // AB -> AA
         truegen_right==atgpk[3]) // AB -> BB
        return(0.5);
      if(truegen_right==atgpk[1]) // AB -> AB
        return(0.0);
      if(truegen_right==atgpk[2]) // AB -> BA
        return(1.0);
    }
    else if(truegen_left == atgpk[2]) {
      if(truegen_right==atgpk[0] || // BA -> AA
         truegen_right==atgpk[3]) // BA -> BB
        return(0.5);
      if(truegen_right==atgpk[2]) // BA -> BA
        return(0.0);
      if(truegen_right==atgpk[1]) // BA -> AB
        return(1.0);
    }
    else if(truegen_left == atgpk[3]) {
      if(truegen_right==atgpk[0]) // BB -> AA
        return(1.0);
      if(truegen_right==atgpk[1] || // BB -> AB
         truegen_right==atgpk[2])   // BB -> BA
        return(0.5);
      if(truegen_right==atgpk[3]) // BB -> BB
        return(0.0);
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double nrecXmale(TrueGenotype truegen_left, TrueGenotype truegen_right)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // AA -> AA
        return(0.0);
      if(truegen_right==atg[3]) // AA -> BB
        return(1.0);
    }
    else if(truegen_left == atg[3]) {
      if(truegen_right==atg[3]) // BB -> BB
        return(0.0);
      if(truegen_right==atg[0]) // BB -> AA
        return(1.0);
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double nrecXfemaleforw(TrueGenotype truegen_left, TrueGenotype truegen_right)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[0]) {
      if(truegen_right==atg[0]) // AA -> AA
        return(0.0);
      if(truegen_right==atg[1]) // AA -> AB
        return(1.0);
    }
    else if(truegen_left == atg[1]) {
      if(truegen_right==atg[1]) // AB -> AB
        return(0.0);
      if(truegen_right==atg[0]) // AB -> AA
        return(1.0);
    }

    throw new Exception("inputs not among the possible true genotypes");
  }

  double nrecXfemalerev(TrueGenotype truegen_left, TrueGenotype truegen_right)
  {
    alias all_true_geno_X atg;

    if(truegen_left==atg[2]) {
      if(truegen_right==atg[2]) // AB -> AB
        return(0.0);
      if(truegen_right==atg[3]) // AB -> BB
        return(1.0);
    }
    else if(truegen_left == atg[3]) {
      if(truegen_right==atg[3]) // BB -> BB
        return(0.0);
      if(truegen_right==atg[2]) // BB -> AB
        return(1.0);
    }

    throw new Exception("inputs not among the possible true genotypes");
  }
}
