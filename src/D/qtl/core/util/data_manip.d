/**
 * Data manipulation
 */

module qtl.core.util.data_manip;

import std.conv, std.stdio, std.string, std.exception, std.path, std.file;
import qtl.core.primitives;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;

import std.typecons;
import std.algorithm;
import std.regex;
import std.variant;


// return boolean vector indicating whether any
bool[] is_any_phenotype_missing(T)(Phenotype!T[][] pheno)
{
  auto ret = new bool[](pheno.length);
  foreach(i; 0..pheno.length) {
    ret[i] = false;
    foreach(j; 0..pheno[i].length) {
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

// omit individuals from genotype matrix
GenotypeCombinator[][] omit_ind_from_genotypes(GenotypeCombinator[][] geno, bool[] to_omit)
{
  if(geno.length != to_omit.length)
    throw new Exception("no. individuals in geno (" ~ to!string(geno.length) ~
                        ") doesn't match length of to_omit (" ~ to!string(to_omit.length) ~ ")");

  GenotypeCombinator[][] ret;

  foreach(i; 0..to_omit.length) {
    if(!to_omit[i])
      ret ~= geno[i];
  }

  return ret;
}

unittest {
  writeln("Unit test " ~ __FILE__);

  // Phenotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];

  // find individuals with missing phenotypes
  auto has_missing = is_any_phenotype_missing(pheno);

  // count them
  auto n_missing = count(has_missing, true);

  // check results
  assert(n_missing == 4);
  auto has_missing_shouldbe = new bool[](120);
  foreach(ref val; has_missing_shouldbe) val=false;
  has_missing_shouldbe[29] = true;
  has_missing_shouldbe[71] = true;
  has_missing_shouldbe[75] = true;
  has_missing_shouldbe[76] = true;
  assert(has_missing == has_missing_shouldbe);

  // omit individuals with missing phenotypes
  auto pheno_rev = omit_ind_from_phenotypes(pheno, has_missing);
  assert(pheno_rev.length == 116);
  auto has_missing_rev = is_any_phenotype_missing(pheno_rev);
  assert(count(has_missing_rev, true) == 0);

  // First read symbol information (the GenotypeCombinators)
  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
  writeln("reading ",symbol_fn);
  auto f = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(f);
  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"regression","test_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto f1 = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(f1, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  assert(genotype_matrix.length == 120);
  // omit individuals with missing phenotypes
  auto genotype_matrix_rev = omit_ind_from_genotypes(genotype_matrix, has_missing);
  assert(genotype_matrix_rev.length == 116);
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

unittest {
  writeln("Test create_phenotype_batches");

  Phenotype!double pheno[][];
  auto pdbl = [ ["0.0", "0.0", "0.0", "0.0", "0.0"],  
                ["0.0", "0.0", "0.0", "0.0", "0.0"], 
                [  "-", "0.0", "0.0",   "-", "0.0"], 
                ["0.0", "0.0", "0.0", "0.0", "0.0"], 
                ["0.0", "0.0", "0.0", "0.0", "0.0"], 
                ["0.0", "0.0", "0.0", "0.0", "0.0"], 
                [  "-",   "-",   "-",   "-",   "-"], 
                ["0.0", "0.0", "0.0", "0.0", "0.0"], 
                [  "-",   "-", "0.0",   "-",   "-"], 
                ["0.0", "0.0", "0.0", "0.0", "0.0"], 
                ["0.0",   "-",   "-", "0.0",   "-"], 
                ["0.0",   "-",   "-", "0.0",   "-"], 
                ["0.0", "0.0", "0.0", "0.0", "0.0"]];

  foreach(line; pdbl) {
      Phenotype!double[] ps = std.array.array(map!((a) {return set_phenotype!double(a);})(line));
      pheno ~= ps;
  }

  // original version
  auto batches = create_phenotype_batches(pheno);
  assert(batches.length == 3);
  assert(batches[0] == [0,3]);
  assert(batches[1] == [1,4]);
  assert(batches[2] == [2]);
  writeln("batches: ", batches);

  // test creation of patterns
  assert(create_missing_phenotype_pattern(pheno, 0) == "2|6|8");
  assert(create_missing_phenotype_pattern(pheno, 1) == "6|8|10|11");
  assert(create_missing_phenotype_pattern(pheno, 2) == "6|10|11");

  // hash-based version
  auto batches2 = create_phenotype_batches_hash(pheno);
  writeln("batches2: ", batches2);
}

