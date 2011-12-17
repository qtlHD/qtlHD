/**
 * Read tabular qtlHD files (.qtab)
 */

module qtl.plugins.qtab.read_qtab;

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
 * Read a tabular qtlHD symbol genotype line
 */

Tuple!(string, string[]) parse_symbol_genotype_qtab(string line) {
  auto fields1 = split(line," ");
  auto fields = std.array.array(map!"strip(a)"(fields1));  // <- note conversion to array
  auto ind = (fields.length > 0 ? strip(fields[0]) : null);
  auto phenotypes = (fields.length > 1 ? fields[1..$] : null);
  return tuple(ind,phenotypes);
}

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
 * Low level qtlHD qtab parser. Read the file in sections, and
 * delegate to the appropriate readers.
 */
void read_qtab(string fn) {
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  string buf;
  while (f.readln(buf)) {
    if (strip(buf) == "# --- Symbol Genotype begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- Symbol Genotype end")
           break;
        auto res = parse_symbol_genotype_qtab(buf);
        auto ind = res[0];
        auto ps  = res[1];
        writeln(ind,"\t",ps);
      }
    }
    if (strip(buf) == "# --- Data Phenotype begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- Data Phenotype end")
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
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data"));
  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
  writeln("reading ",symbol_fn);
  read_qtab(symbol_fn);
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  read_qtab(pheno_fn);
}

