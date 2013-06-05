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

import std.exception, std.path, std.file;
import qtl.core.primitives;
import qtl.core.phenotype;

import std.typecons;
import std.algorithm : map;


Phenotype!T set_phenotype(T)(in string s) {
  // writeln(s);
  Phenotype!T p;
  if(s == "NA" || s == "-"){
    p.value = to!T(PHENOTYPE_NA);
  }else{
    if(s.indexOf(".") != -1){  // FIXME: this should only be for floats
      p.value = to!T(s);
    }else{
      p.value = to!T(s~".0");
    }
  }
  return p;
}

/**
 * Check whether a phenotype is missing
 */
bool isNA(T)(Phenotype!T phe) { 
  return(phe.value == PHENOTYPE_NA);
}

// return boolean vector of size individuals indicating whether a 
// phenotype is missing (true)
bool[] individuals_missing_a_phenotype(T)(Phenotype!T[][] pheno)
{
  auto ret = new bool[pheno.length];
  ret[] = false;
  foreach(i; 0..pheno.length) {
    // walk all individuals
    foreach(j, p_list; pheno[i]) {
      // walk all phenotypes
      if(isNA(pheno[i][j])) ret[i] = true;
      break;
    }
  }
  return ret;
}

// omit individuals from phenotype data
Phenotype!T[][] omit_ind_from_phenotypes(T)(Phenotype!T[][] pheno, bool[] to_omit)
{
  if(pheno.length != to_omit.length)
    throw new Exception("no. individuals in pheno (" ~ to!string(pheno.length) ~
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  Phenotype!T[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= pheno[i];
  }

  return ret;
}

// batches of phenotypes with common missing data patterns
size_t[][] create_phenotype_batches(T)(Phenotype!T[][] pheno)
{
  size_t[][] phenotype_batches = [ [0] ];

  size_t[] first_pattern = [];
  size_t[][] patterns;

  if(pheno[0].length == 1) return phenotype_batches;

  foreach(i; 0..pheno.length) {
    if(isNA(pheno[i][0]))
      first_pattern ~= i;
  }
  patterns ~= first_pattern;

  foreach(j; 1..pheno[0].length) {
    size_t[] this_pattern = [];

    foreach(i; 0..pheno.length) {
      if(isNA(pheno[i][j]))
        this_pattern ~= i;
    }

    bool found = false;
    foreach(i; 0..patterns.length) {
      if(this_pattern == patterns[i]) {
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
string create_missing_phenotype_pattern(T)(Phenotype!T[][] pheno, size_t pheno_column)
{
  string pattern;

  foreach(i; 0..pheno.length) {
    if(isNA(pheno[i][pheno_column])) {
      if(pattern == "") pattern = to!string(i);
      else pattern ~= "|" ~ to!string(i);
    }
  }

  return pattern;
}


// batches of phenotypes with common missing data patterns
size_t[][string] create_phenotype_batches_hash(T)(Phenotype!T[][] pheno)
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


