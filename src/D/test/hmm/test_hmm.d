/**
 * test reading data and then apply HMM
 */
 
module test.hmm.test_hmm;

import std.math, std.stdio, std.path;
import std.exception, std.conv;

import qtl.core.primitives;
import qtl.core.phenotype, qtl.core.chromosome, qtl.core.genotype;
import qtl.plugins.qtab.read_qtab;

unittest {
  writeln("------------------");
  writeln("unit test test_hmm");
  writeln("------------------");

  // Symbol and genotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","test","data", "input", "listeria_qtab"));
  auto symbol_fn = to!string(buildPath(dir,"listeria_symbol.qtab"));
  // First read symbol information (the GenotypeCombinators)
  writeln("reading ",symbol_fn);
  auto fs = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(fs);

  // Test working of symbols
  assert(symbols.decode("A") == symbols.decode("AA"));
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("H")) == "[(0,1), (1,0)]");
  //assert(to!string(symbols.decode("HorA")) == "[(0,0), (0,1)]", to!string(symbols.decode("HorA")));

  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"listeria_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto fg = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(fg, symbols);
  auto individuals = ret[0];
  auto genotype_matrix = ret[1];

  // Show the first individual and genotypes
  // writeln(individuals.list[0].name,genotype_matrix[0]);
  // by symbol
  assert(genotype_matrix[0][0] == symbols.decode("B"));
  assert(genotype_matrix[0][3] == symbols.decode("H"));
  // by founders
  assert(genotype_matrix[0][0].list[0].homozygous == true);
  assert(genotype_matrix[0][0].list[0].founders[0] == 1);
  assert(genotype_matrix[0][0].list[0].founders[1] == 1);
  assert(genotype_matrix[0][3].list[0].heterozygous == true);
  assert(genotype_matrix[0][3].list[0].founders[0] == 0);
  assert(genotype_matrix[0][3].list[0].founders[1] == 1);
  // assert(genotype_matrix[0][3].list[0] == 1);

  // reading phenotypes
  auto pheno_fn = to!string(buildPath(dir,"listeria_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];
  // 1st ind, 1st phenotype
  assert(pheno[0][0].value == 118.317);
  // 3rd ind, 1st phenotype
  assert(pheno[2][0].value == 194.917);

  // Marker map reader
  auto marker_map_fn = to!string(buildPath(dir,"listeria_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);
  assert(markers[0].name == "D10M44");
  assert(markers[0].chromosome.name == "1");
  assert(markers[3].position == 40.4136);

}



