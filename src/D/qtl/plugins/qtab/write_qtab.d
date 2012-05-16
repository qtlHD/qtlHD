/**
 * Write tabular qtlHD files (.qtab)
 *
 * See ./doc/input/qtab.md for a definition of the qtab standard
 */

module qtl.plugins.qtab.write_qtab;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.individual;
import qtl.plugins.csv.read_csv;  // for testing
import example.genotype_examples; // for testing

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

immutable string VER = "0.1";
string ID = "qtlHD-in-" ~ VER;

/**
 * Write a qtab (tabular) qtlHD symbol section
 */

void write_symbol_qtab(File f, string descr, ObservedGenotypes observed) {
  f.writeln("# --- ",ID," Symbol ",descr);
  f.writeln("# --- Genotype Symbol begin");
  foreach(symbol ; observed.list) {
    f.write(symbol.toEncodingString," as ");
    f.write(symbol.toTrueGenotypes);
    f.writeln();
  }
  f.writeln("# --- Genotype Symbol end");
}

/**
 * Write a qtab (tabular) qtlHD marker map section
 */

void write_marker_map_qtab(File f, string descr, Marker[] markers) {
  f.writeln("# --- ",ID," Location ",descr);
  f.writeln("# --- Data Location begin");
  f.writeln("#\tChr\tPos");
  foreach(m ; markers.list) {
    f.write(m.name);
    f.write("\t",m.chromosome.name);
    f.write("\t",m.position);
    f.writeln();
  }
  f.writeln("# --- Data Location end");
}


/**
 * Write a tabular qtlHD genotype section
 */

void write_genotype_qtab(G)(File f, string descr, in Individuals individuals, in Marker[] markers, in G[][] genotypes) {
  f.writeln("# --- qtlHD-in-0.1 Genotype " ~ descr);
  f.writeln("# --- Data Genotype begin");
  // write markers as a comment
  f.write("#");
  foreach(m; markers) {
    f.write("\t",m.name);
  }
  f.writeln;
  foreach(i, ind; individuals.list) {
    f.write(ind.name);
    foreach(g ; genotypes[i]) {
      f.write("\t",g.toEncoding);
    }
    f.writeln;
  }
  f.writeln("# --- Data Genotype end");
}

/**
 * Write a tabular qtlHD phenotype section
 */

void write_phenotype_qtab(P)(File f, string descr, in Individuals individuals, in string[] phenotypenames, in P[][] phenotypes) {
  f.writeln("# --- qtlHD-in-0.1 Phenotype " ~ descr);
  // first we write the names and types
  f.writeln("# --- Type Phenotype begin");
  foreach(n; phenotypenames) {
    f.writeln(n,"\tFloat");
  }
  f.writeln("# --- Type Phenotype end");

  f.writeln("# --- Data Phenotype begin");
  // write phenotype names as a comment
  f.write("#");
  foreach(n; phenotypenames) {
    f.write("\t",n);
  }
  f.writeln;
  foreach(i, ind; individuals.list) {
    f.write(ind.name);
    foreach(p; phenotypes[i]) {
      f.write("\t",p.toString);
    }
    f.writeln;
  }
  f.writeln("# --- Data Phenotype end");
}

unittest {
  // in this unit test we read a CSV file, and write it into 
  // valid qtab formatted files in the regression directory. Note the read_qtab unit tests
  // read these same files, for testing
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));
  auto fn = to!string(buildPath(dir,"input","listeria.csv"));
  auto data = new ReadSimpleCSV!(F2,ObservedF2)(fn);
  // write the phenotype file
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab.new"));
  auto f = File(pheno_fn,"w");
  write_phenotype_qtab(f, "Test", data.individuals, data.phenotypenames, data.phenotypes);
  f.close();
  writeln("Wrote ", pheno_fn);
  // write the marker map file 
  auto marker_map_fn = to!string(buildPath(dir,"regression","test_marker_map.qtab.new"));
  f = File(marker_map_fn,"w");
  write_marker_map_qtab(f, "Test", data.markers);
  f.close();
  writeln("Wrote ", marker_map_fn);
  // write the genotype file 
  auto geno_fn = to!string(buildPath(dir,"regression","test_genotype.qtab.new"));
  f = File(geno_fn,"w");
  write_genotype_qtab(f, "Test", data.individuals, data.markers, data.genotypes);
  f.close();
  writeln("Wrote ", geno_fn);
  // write the symbol file 
  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab.new"));
  f = File(symbol_fn,"w");
  write_symbol_qtab(f, "Test", data.symbols.symbols);
  f.close();
  writeln("Wrote ", symbol_fn);
}
