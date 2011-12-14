/**
 * Write tabular qtlHD files (.qtab)
 */

module qtl.plugins.output.write_qtab;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.plugins.input.read_csv;  // for testing

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

/**
 * Write a tabular qtlHD genotype file. 
 */

void write_genotype_qtab(File f, string descr, in Individuals individuals, in Marker[] markers) {
  f.writeln("# --- qtlHD-in-0.1 Genotype " ~ descr);
  f.writeln("# --- Data Genotypes begin");
  // write markers as a comment
  f.write("#");
  foreach(m; markers) {
    f.write("\t",m.name);
  }
  f.writeln;
  foreach(i, ind; individuals.list) {
    f.write(ind.name);
    // foreach(p; phenotypes[i]) {
    //  f.write("\t",p.toString);
    // }
    f.writeln;
  }
  f.writeln("# --- Data Genotypes end");
}

/**
 * Write a tabular qtlHD phenotype file. 
 */

void write_phenotype_qtab(P)(File f, string descr, in Individuals individuals, in string[] phenotypenames, in P[][] phenotypes) {
  f.writeln("# --- qtlHD-in-0.1 Phenotype " ~ descr);
  // first we write the names and types
  f.writeln("# --- Type Phenotypes begin");
  foreach(n; phenotypenames) {
    f.writeln(n,"\tFloat");
  }
  f.writeln("# --- Type Phenotypes end");

  f.writeln("# --- Data Phenotypes begin");
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
  f.writeln("# --- Data Phenotypes end");
}

unittest {
  // in this unit test we read a CSV file, and write it into 
  // valid qtab formatted files. Note the read_qtab unit tests
  // read these same files in for testing
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data"));
  auto fn = to!string(buildPath(dir,"input","listeria.csv"));
  auto data = new ReadSimpleCSV!(F2,ObservedF2)(fn);
  // write the phenotype file
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab"));
  auto f = File(pheno_fn,"w");
  write_phenotype_qtab(f, "Test", data.individuals, data.phenotypenames, data.phenotypes);
  f.close();
  // write the genotype file
  auto geno_fn = to!string(buildPath(dir,"regression","test_genotype.qtab"));
  f = File(geno_fn,"w");
  write_genotype_qtab(f, "Test", data.individuals, data.markers);
  f.close();
}
