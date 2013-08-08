/**
 * Phenotype module
 *
 * PhenotypeMatrix holds phenotypes (cols) against individuals (rows)
 */

module qtl.core.phenotype;

import std.conv;
import std.stdio;
import std.string;
import std.array;
import std.typecons;
import std.algorithm;
import std.exception;
import std.path;
import std.file;
import std.math;

import qtl.core.primitives;
import qtl.core.phenotype;
import qtl.core.data.matrix;

immutable PHENOTYPE_NA = double.nan; 

/**
 * AnyPhenotype is the most primitive representation of a phenotype. The type
 * can be any type T (normally a double, but can potentially be any Object).
 *
 * Note the primitive should be small as small as possible, there may be many
 * phenotypes! Therefore it is a struct.
 *
 * Note we do not use this facility other than outputting 'NA' when missing.
 */

struct AnyPhenotype(T) {
  T value;
  
  /// String representation of phenotype.
  const string toString(){
    if(to!double(value) != PHENOTYPE_NA){
      return to!string(value);
    }else{
      return "NA";
    }
  }
}

alias AnyPhenotype!double Phenotype ;

/**
 * PhenotypeMatrix holds phenotypes (cols) against individuals (rows)
 */

alias Phenotype[][] PhenotypeMatrix; // = new double[][](n_ind,n_phe);

Phenotype set_phenotype(in string s) {
  // writeln(s);
  Phenotype p;
  if(s == "NA" || s == "-"){
    p.value = PHENOTYPE_NA;
  }else{
    if(s.countUntil(".") != -1){  // FIXME: this is only for float/double
      p.value = to!double(s);
    }else{
      p.value = to!double(s~".0");
    }
  }
  return p;
}

/**
 * Check whether a phenotype is missing
 */
bool isNA(Phenotype p) { 
  return(isNaN(p.value));
}

// Comparison, including missingness
bool isSame(Phenotype p1, Phenotype p2) {
  if(isNaN(p1.value) && isNaN(p2.value)) return true;
  if(!isNaN(p1.value) && !isNaN(p2.value) && p1.value == p2.value) return true;
  return false;
}

// return boolean vector of size individuals indicating whether a 
// phenotype is missing (true)
bool[] individuals_missing_a_phenotype(Phenotype[][] phenotype_matrix)
{
  return filter_matrix_by_row_2bool!Phenotype(phenotype_matrix, (p) => isNA(p));
}

// omit individuals from phenotype data
Phenotype[][] omit_ind_from_phenotypes(Phenotype[][] pheno, bool[] to_omit)
{
  if(pheno.length != to_omit.length)
    throw new Exception("no. individuals in pheno (" ~ to!string(pheno.length) ~
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  Phenotype[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= pheno[i];
  }

  return ret;
}

// select a subset of phenotype columns
Phenotype[][] subset_phenotype_columns(Phenotype[][] pheno, size_t[] pheno_columns)
{
  foreach(col; pheno_columns)
    if(col < 0 || col >= pheno[0].length)
      throw new Exception("pheno_columns outside of allowable range");

  Phenotype[][] ret;

  foreach(pherow; pheno) {
    Phenotype[] thisrow;
    foreach(col; pheno_columns)
      thisrow ~= pherow[col];
    ret ~= thisrow;
  }

  return ret;
}

// batches of phenotypes with common missing data patterns
size_t[][] create_phenotype_batches(Phenotype[][] pheno)
{
  size_t[][] phenotype_batches = [ [0] ];

  size_t[] first_pattern = [];
  size_t[][] patterns;

  if(pheno[0].length == 1) return phenotype_batches;

  foreach(i, p; pheno) {
    if(isNA(p[0]))
      first_pattern ~= i;
  }
  patterns ~= first_pattern;

  foreach(j; 1..pheno[0].length) {
    size_t[] this_pattern = [];

    foreach(i, p; pheno) {
      if(isNA(p[j]))
        this_pattern ~= i;
    }

    bool found = false;
    foreach(i, pat; patterns) {
      if(this_pattern == pat) {
        found = true;
        phenotype_batches[i] ~= j;
      }
    }
    if(!found) {
      phenotype_batches ~= [j];
      patterns ~= this_pattern;
    }
  }

  return phenotype_batches;
}

// create string with pattern of missing values for a phenotype column
//    if individuals 2, 6, 8 are missing, the output will be "2|6|8"
string create_missing_phenotype_pattern(Phenotype[][] pheno, size_t pheno_column)
{
  string pattern;

  foreach(i, p; pheno) {
    if(isNA(p[pheno_column])) {
      if(pattern == "") pattern = to!string(i);
      else pattern ~= "|" ~ to!string(i);
    }
  }

  return pattern;
}


// batches of phenotypes with common missing data patterns
size_t[][string] create_phenotype_batches_hash(Phenotype[][] pheno)
{
  auto pattern = create_missing_phenotype_pattern(pheno, 0);
  auto phenotype_batches = [pattern: [cast(size_t)0]];
  if(pheno[0].length == 1) return phenotype_batches;

  foreach(i; 1..pheno[0].length) {
    pattern = create_missing_phenotype_pattern(pheno, i);
    phenotype_batches[pattern] ~= i;
  }

  return phenotype_batches;
}

// ---- Phenotype vector to double[]
T[] get_values(T)(AnyPhenotype!T[] input)
{
  return std.array.array(map!"a.value"(input));
}

// ---- Phenotype matrix to double[][]
T[][] get_values(T)(AnyPhenotype!T[][] input)
{
  return std.array.array(map!(a => get_values(a))(input));
}
