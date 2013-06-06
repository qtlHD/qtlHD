module test.qtl.core.phenotype;

import std.conv;
import std.stdio;
import std.string;

import std.exception, std.path, std.file;
import qtl.core.primitives;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;

import std.typecons;
import std.algorithm : map;


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
  auto has_missing = individuals_missing_a_phenotype(pheno);

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
  auto has_missing_rev = individuals_missing_a_phenotype(pheno_rev);
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
  auto batches_hash = create_phenotype_batches_hash(pheno);
  writeln("batches_hash: ", batches_hash);
  assert(batches_hash.length == 3);
  assert(batches_hash["2|6|8"] == [0,3]);
  assert(batches_hash["6|8|10|11"] == [1,4]);
  assert(batches_hash["6|10|11"] == [2]);
}

