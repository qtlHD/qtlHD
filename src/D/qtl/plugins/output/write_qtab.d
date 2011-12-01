/**
 * Write tabular qtlHD files (.qtab)
 */

module qtl.plugins.output.write_qtab;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

/**
 * Write a tabular qtlHD phenotype file. 
 */

void write_phenotype_qtab(P)(in File f, in Individuals individuals, in P[][] phenotypes) {
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  // auto fn = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data","input","listeria.csv"));
  auto f = File("test_phenotype.qtab","w");
  write_phenotype_qtab!double(f, null, null);
  
  f.close();
}
