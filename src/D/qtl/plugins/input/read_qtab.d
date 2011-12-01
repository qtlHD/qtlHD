/**
 * Read tabular qtlHD files (.qtab)
 */

module qtl.plugins.input.read_qtab;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;
import std.typecons;
import std.algorithm;

/**
 * Read a tabular qtlHD phenotype line
 */

Tuple!(string, string[]) parse_phenotype_qtab(string line) {
  auto fields1 = split(line,"\t");
  auto fields = std.array.array(map!"strip(a)"(fields1));  // <- note conversion to array
  auto ind = (fields.length > 0 ? strip(fields[0]) : null);
  auto phenotypes = (fields.length > 1 ? fields[1..$] : null);
  return tuple(ind,phenotypes);
}

/**
 * Low level qtlHD qtab parser
 */
void read_qtab(string fn) {
  auto f = File(fn,"r");
  scope(exit) f.close();
  string buf;
  while (f.readln(buf)) {
    if (strip(buf) == "# --- Data Phenotypes begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- Data Phenotypes end")
           break;
        auto res = parse_phenotype_qtab(buf);
        auto ind = res[0];
        auto ps  = res[1];
        writeln(ind,"\t",ps);
      }
    }
  }
}

unittest {
  writeln("Unit test " ~ __FILE__);
  read_qtab("test_phenotype.qtab");

}

