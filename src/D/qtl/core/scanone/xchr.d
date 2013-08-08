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
    cross_dir = new int[][](sex.length, 1);
    foreach(i, val; tmp)
      cross_dir[i][0] = cast(int)(val.value);
  }
  else {
    throw new Exception("cross type " ~ cross.cross_type ~ " not supported.");
  }

  return tuple(is_female, cross_dir);
}

// ---- Test if all male
bool all_male(bool[] is_female)
{
  foreach(val; is_female) {
    if(val) return(false);
  }
  return(true); // if is_female.length == 0, still return true
}

// ---- Test if all female
bool all_female(bool[] is_female)
{
  foreach(val; is_female) {
    if(!val) return(false);
  }
  return(true); // if is_female.length == 0, still return true
}

// ---- Test if all one sex
bool all_same_sex(bool[] is_female)
{
  foreach(val; is_female) {
    if(val != is_female[0]) return(false);
  }
  return(true); // if is_female.length == 0, still return true
}

// ---- Test if all same cross direction
bool all_same_cross_direction(int [][] cross_direction)
{
  foreach(val; cross_direction) {
    if(val[0] != cross_direction[0][0]) return(false);
  }
  return(true);
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
  if(all_same_sex(is_female)) { // all females or all males, no covariate
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

  auto all_same_sex = all_same_sex(is_female);
  auto all_same_cross_direction = all_same_cross_direction(cross_direction);

  if(all_same_sex) { // all same sex
    if(all_same_cross_direction || !is_female[0]) { // all males, or all females of same cross direction
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