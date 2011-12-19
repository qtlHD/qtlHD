/**
 * Read tabular qtlHD files (.qtab)
 */

module qtl.plugins.qtab.read_qtab;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;
import std.typecons;
import std.algorithm;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

/**
 * Read a tabular qtlHD symbol genotype line - splits the line on
 * spaces, followed by looking for 'as'. E.g.
 *
 *   AC CA as 0,2 2,0
 *
 * gets returned as
 *
 *   Tuple(["AC","CA"],["0,2","2,0"])
 */

Tuple!(string[], string[]) parse_symbol_genotype_qtab(string line) {
  auto fields1 = split(line," ");
  auto fields = std.array.array(map!"strip(a)"(fields1));  // <- note conversion to array
  auto res = find(fields,"as");
  auto genotypes = (res.length > 1 ? res[1..$] : null);
  auto symbols = fields[0..$-genotypes.length-1];
  return tuple(symbols,genotypes);
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
 * delegate to the appropriate readers. Values are allocated
 * in memory and returned in a Tuple.
 */
Tuple!(P[][]) read_qtab(P)(string fn) {
  P ret_phenotypes[][];  // return matrix
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  string buf;
  while (f.readln(buf)) {
    if (strip(buf) == "# --- Symbol Genotype begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- Symbol Genotype end")
           break;
        auto res = parse_symbol_genotype_qtab(buf);
        auto symbols = res[0];
        auto genotypes = res[1];
        writeln(symbols,"\t",genotypes);
      }
    }
    if (strip(buf) == "# --- Data Phenotype begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- Data Phenotype end")
           break;
        if (buf[0] == '#') continue;
        auto res = parse_phenotype_qtab(buf);
        auto ind = res[0];
        auto fields  = res[1];
        writeln(ind,"\t",fields);
        P[] ps = std.array.array(map!((a) {return set_phenotype!double(a);})(fields));
        ret_phenotypes ~= ps;
      }
    }
  }
  return tuple(ret_phenotypes);
}

unittest {
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..","..","test","data"));
  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
  writeln("reading ",symbol_fn);
  read_qtab!(Phenotype!double)(symbol_fn);
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];
  // 1st ind, 1st phenotype
  assert(pheno[0][0].value == 118.317);
  // 3rd ind, 1st phenotype
  assert(pheno[2][0].value == 194.917);
}

