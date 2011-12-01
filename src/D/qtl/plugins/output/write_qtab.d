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
 * Write a tabular qtlHD phenotype file. 
 */

void write_phenotype_qtab(P)(File f, in Individuals individuals, in string[] phenotypenames, in P[][] phenotypes) {
  // first we write the names and types
  f.writeln("# --- Type Phenotypes begin");
  foreach(n; phenotypenames) {
    f.writeln(n,"\tFloat");
  }
  f.writeln("# --- Type Phenotypes end");

  f.writeln("# --- Data Phenotypes begin");
  // write markers as a comment
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
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto fn = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data","input","listeria.csv"));
  auto data = new ReadSimpleCSV!(F2,ObservedF2)(fn);
  auto f = File("test_phenotype.qtab","w");
  write_phenotype_qtab(f, data.individuals, data.phenotypenames, data.phenotypes);
  f.close();
}
