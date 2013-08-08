/**
 * X chromosome related functions
 */

module qtl.core.scanone.xchr;

import std.algorithm;
import std.stdio;
import std.string;
import std.conv;
import std.typecons;
import std.exception;

import qtl.core.chromosome;
import qtl.core.genotype;
import qtl.core.individual;
import qtl.core.marker;
import qtl.core.phenotype;
import qtl.core.primitives;
import qtl.core.hmm.cross;
import qtl.core.data.matrix;

// ---- Get sex and cross direction from phenotypes
Tuple!(bool [], int [][]) get_sex_and_cross(Cross cross, string sexcolumn, string crosscolumn, 
                                            string[] phenames, Phenotype[][] pheno_matrix)
{
  auto sex = get_phenotype(sexcolumn, phenames, pheno_matrix);
  auto is_female = new bool[](sex.length);
  foreach(i, val; sex)
    is_female[i] = (val.value == 0.0);

  int[][] cross_dir;

  if(cross.cross_type == "BC") {
    cross_dir = new int[][](sex.length, 0);
  }
  else if(cross.cross_type == "F2") {
    auto tmp = get_phenotype(crosscolumn, phenames, pheno_matrix);
    foreach(i, val; tmp)
      cross_dir[i][0] = cast(int)(val.value);
  }
  else {
    throw new Exception("cross type " ~ cross.cross_type ~ " not supported.");
  }

  return tuple(is_female, cross_dir);
}


// ---- Special covariates for QTL analysis of the X chromosome
double[][] xchr_covar(Cross cross, bool[] is_female, int[][] cross_direction)
{
  if(cross.cross_type == "BC") {
    return xchr_covar_bc(is_female);
  }
  else if(cross.cross_type == "F2") {
    return xchr_covar_f2(is_female, cross_direction);
  }
  else {
    throw new Exception("cross type " ~ cross.cross_type ~ " not supported.");
  }
}

// ---- Version that takes additive covariate matrix as input and gives combined matrix as output
double [][] xchr_covar(Cross cross, bool[] is_female, int[][] cross_direction, double [][] addcovar)
{
  auto newcovar = xchr_covar(cross, is_female, cross_direction);

  return(bind_columns(addcovar, newcovar));
}



// ---- Backcross: special covariates for X chr
double [][] xchr_covar_bc(bool[] is_female)
{
  if(is_female.length < 2) { // 0 or 1 individual
    auto ret = new double[][](is_female.length, 0);
    return ret;
  }

  auto all_same = true;
  foreach(i; 1..is_female.length) {
    if(is_female[i] != is_female[0]) {
      all_same = false;
      break;
    }
  }
  
  if(all_same) { // all females or all males, no covariate
    auto ret = new double[][](is_female.length, 0);
    return ret;
  }

  // some males and some females
  auto x_covar = new double[][](is_female.length, 1);
  foreach(i, isF; is_female)
    x_covar[i][0] = cast(double)(!isF); // female = 0, male = 1

  return x_covar;
}

// ---- Intercross: special covariates for X chr
double [][] xchr_covar_f2(bool[] is_female, int[][] cross_direction)
{
  if(is_female.length != cross_direction.length)
    throw new Exception("is_female and cross_direction have different lengths: "
                        ~ to!string(is_female.length) ~ " " ~ to!string(cross_direction.length));

  auto n_ind = is_female.length;

  if(n_ind < 2) { // 0 or 1 individual
    auto ret = new double[][](n_ind, 0);
    return ret;
  }

  auto all_same_sex = true;
  foreach(i; 1..n_ind) {
    if(is_female[i] != is_female[0]) {
      all_same_sex = false;
      break;
    }
  }
  
  auto all_same_cross_direction = true;
  foreach(i; 1..n_ind) {
    if(cross_direction[i][0] != cross_direction[0][0]) {
      all_same_cross_direction = false;
      break;
    }
  }
  
  if(all_same_sex) { // all same sex
    if(!is_female[0] || all_same_cross_direction) { // all males, or all females of same cross direction
      auto ret = new double[][](n_ind, 0);
      return ret;
    }
    else { // all females but both directions: covar = (sex)
      auto ret = new double[][](n_ind, 1);
      foreach(i, val; cross_direction)
        ret[i][0] = cast(double)val[0];
      return ret;
    }
  }
  else { // both sexes
    if(all_same_cross_direction) { // both sexes, one direction: covar = (sex)
      auto ret = new double[][](n_ind, 1);
      foreach(i, isF; is_female)
        ret[i][0] = cast(double)(!isF);
      return ret;
    }
    else { // both sexes and both crosses: covar = (male/female, cross direction for females)
      auto ret = new double[][](n_ind, 2);
      foreach(i, val; cross_direction) {
        if(is_female[i]) {
          ret[i] = [0.0, cast(double)val[0]];
        }
        else {
          ret[i] = [1.0, 0.0];
        }
      }
      return ret;
    }
  }
}

// ---- revise layout of QTL genotype probabilities for X chromosome
//double[][][] revise_X_genoprobs(Cross cross, int[] sex, int[] crossdir, double [][][] probs)
//{
//}