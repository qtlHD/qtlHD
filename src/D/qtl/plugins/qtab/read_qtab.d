/**
 * Read tabular qtlHD files (.qtab)
 *
 * See ./doc/input/qtab.md for a definition of the qtab standard
 */

module qtl.plugins.qtab.read_qtab;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;
import std.typecons;
import std.algorithm;
import std.regex;
alias std.string.split str_split;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.individual;

/**
 * Low level parsers, str_splits a tab delimited line into fields,
 * and removes everything after a remark '#' symbol. Optionally
 * strips each field of white space.
 */


string[] str_split_line(string line,string split_with="\t") {
  // strip remark and leading whitespace
  auto res = splitter(chomp(line), regex("\\s+#\\s"));
  return str_split(to!string(std.array.array(res)[0]),split_with);
}

string[] str_split_line_on_spaces(string line) {
  return str_split_line(line," ");
}

string[] str_split_line_strip(string line) {
  auto fields1 = str_split_line(line);
  return std.array.array(map!"strip(a)"(fields1));
}

unittest {
  writeln("Unit test " ~ __FILE__);
  assert(str_split_line("test  # remark") == ["test"],to!string(str_split_line("test  # remark")));
  assert(str_split_line("test\n") == ["test"],to!string(str_split_line("test")));
  assert(str_split_line("test\ttest") == ["test","test"]);
  assert(str_split_line("test\t test \t ") == ["test"," test "," "]);
  assert(str_split_line_strip("test\t test\t") == ["test","test",""]);
  assert(str_split_line_strip("test\t test\t") == ["test","test",""]);
  assert(str_split_line_strip("test\t test\t# remark") == ["test","test"]);
  assert(str_split_line_strip("test\t test\t # remark") == ["test","test"]);
  assert(str_split_line_on_spaces("test test\t # remark") == ["test","test"]);
  writeln(str_split_line_on_spaces("test test  test \n"));
  assert(str_split_line_on_spaces("test test  test \n") == ["test","test","","test",""]);
}

/** 
 * Parse a line for key values
 */

Tuple!(string, string) parse_line_key_values(string line) {
  auto fields = str_split_line_strip(line);
  auto key = (fields.length > 0 ? fields[0] : null);
  auto value = (fields.length > 1 ? fields[1] : null);
  return tuple(key,value);
}

/**
 * Call function 'call_line' for each line in a QTAB file section
 */
void each_line_in_section(File f, string tag, void delegate (string) call_line) {
  string buf;
  f.rewind();
  while (f.readln(buf)) {
    if (strip(buf) == "# --- "~tag~" begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- "~tag~" end")
           break;
        if (buf[0] == '#') continue;
        writeln("Line: ",buf); 
        call_line(buf);
      }
    }
  }
}

/**
 * Turn a QTAB file section into key-value pairs
 */

string[string] get_section_key_values(File f, string tag) {
  string[string] ret;
  each_line_in_section(f,tag, (line) {
    auto res = parse_line_key_values(line);
    ret[res[0]] = res[1];
  });
  return ret;
}

/**
 * Parse a Set Founder section, and return the results in a Hash
 */

string[string] read_founder_settings_qtab(string fn) {
  
  string[string] ret;
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  return get_section_key_values(f,"Set Founder");
}

unittest {
  // Founder reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));

  auto founder_fn = to!string(buildPath(dir,"input","listeria_qtab","listeria_founder.qtab"));
  writeln("reading ",founder_fn);
  auto info = read_founder_settings_qtab(founder_fn);
  writeln(info);
  assert(info["Cross"] == "F2");
}

/**
 * Read a tabular qtlHD phenotype line
 */

Tuple!(string, string[]) parse_phenotype_qtab(string line) {
  auto fields1 = str_split(line,"\t");
  auto fields = std.array.array(map!"strip(a)"(fields1));  // <- note conversion to array
  auto ind = (fields.length > 0 ? strip(fields[0]) : null);
  auto phenotypes = (fields.length > 1 ? fields[1..$] : null);
  return tuple(ind,phenotypes);
}

Tuple!(string, string, double) parse_marker_qtab(string line) {
  auto fields1 = str_split(line,"\t");
  auto fields = std.array.array(map!"strip(a)"(fields1));  // <- note conversion to array
  auto name = (fields.length > 0 ? strip(fields[0]) : null);
  auto chromosome = (fields.length > 1 ? fields[1] : null);
  auto position = (fields.length > 2 ? to!double(strip(fields[2])) : MARKER_POSITION_UNKNOWN);
  return tuple(name,chromosome,position);
}

/**
 * Read a tabular qtlHD genotype line
 */

Tuple!(string, string[]) parse_genotype_qtab(string line) {
  auto fields1 = str_split(line,"\t");
  auto fields = std.array.array(map!"strip(a)"(fields1));  // <- note conversion to array
  auto ind = (fields.length > 0 ? strip(fields[0]) : null);
  auto data = (fields.length > 1 ? fields[1..$] : null);
  return tuple(ind,data);
}


auto read_marker_map_qtab(Ms)(string fn) {  // Ms is Marker[] (vs Markers)
  Ms ret_ms[];
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  string buf;
  uint id=0;
  while (f.readln(buf)) {
    if (strip(buf) == "# --- Data Location begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- Data Location end")
           break;
        if (buf[0] == '#') continue;
        auto res = parse_marker_qtab(buf);
        auto name = res[0];
        auto cname = res[1];
        auto c = new Chromosome(cname);
        auto pos = res[2];
        auto marker = new Marker(c,pos,name,id);
	id++;
	ret_ms ~= marker;
      }
    }
  }
  return ret_ms;
}

/**
 * Low level qtlHD qtab parser. Read the file in sections, and
 * delegate to the appropriate readers. Values are allocated
 * in memory and returned in a Tuple.
 */
Tuple!(P[][]) read_phenotype_qtab(P)(string fn) {
  P ret_phenotypes[][];  // return matrix
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  string buf;
  while (f.readln(buf)) {
    if (strip(buf) == "# --- ??Phenotype Symbol begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- ??Phenotype Symbol end")
           break;
        auto res = parse_symbol_genotype_qtab(buf);
        auto symbols = res[0];
        auto genotypes = res[1];
        // writeln(symbols,"\t",genotypes);
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
        // writeln(ind,"\t",fields);
        P[] ps = std.array.array(map!((a) {return set_phenotype!double(a);})(fields));
        ret_phenotypes ~= ps;
      }
    }
  }
  return tuple(ret_phenotypes);
}

unittest {
  // Phenotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype!double)(pheno_fn);
  Phenotype!double[][] pheno = p_res[0];
  // 1st ind, 1st phenotype
  assert(pheno[0][0].value == 118.317);
  // 3rd ind, 1st phenotype
  assert(pheno[2][0].value == 194.917);
}

unittest {
  // Marker map reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));
  auto marker_map_fn = to!string(buildPath(dir,"regression","test_marker_map.qtab"));
  writeln("reading ",marker_map_fn);
  auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);
  assert(markers[0].name == "D10M44");
  assert(markers[0].chromosome.name == "1");
  assert(markers[3].position == 40.4136);
}

/**
 * Read a tabular qtlHD Genotype Symbol line - str_splits the line on
 * spaces, followed by looking for 'as'. E.g.
 *
 *   AC CA as 0,2 2,0
 *
 * gets returned as
 *
 *   Tuple(["AC","CA"],["0,2","2,0"])
 */

Tuple!(string[], string[]) parse_symbol_genotype_qtab(string line) {
  auto fields = str_split_line_on_spaces(line);
  writeln(fields);
  auto res = find(fields,"as");
  auto genotypes = (res.length > 1 ? res[1..$] : null);
  auto symbols = fields[0..$-genotypes.length-1];
  return tuple(symbols,genotypes);
}

/**
 * Parse genotype symbol table
 *
 * Returns an ObservedGenotypes container
 */

ObservedGenotypes read_genotype_symbol_qtab(File f, bool phase_known = true) {
  auto observed = new ObservedGenotypes();
  each_line_in_section(f,"Genotype Symbol", 
    (line) {
      auto res = parse_symbol_genotype_qtab(line);
      auto symbol_names = res[0];
      auto genotype_strs = res[1];
      writeln("Symbol: ",symbol_names,"\t",genotype_strs);
      auto combinator = new GenotypeCombinator(symbol_names[0], null, phase_known);
      foreach (s ; symbol_names[1..$]) {
        combinator.add_encoding(s);
      }
      foreach (g ; genotype_strs) { 
        if (g == "None") continue;
        combinator ~= new TrueGenotype(g);
      }
      // writeln(combinator.toEncodings);
      observed ~= combinator;
    });
  // writeln(observed);
  return observed;
}

class SymbolSettings {
}

/**
 * Read set symbol section in qtab file and convert to a SetSymbol
 * object
 */

void do_parse(File f,string tag) {
  f.rewind();
}

SymbolSettings read_set_symbol_qtab(File f) {
  do_parse(f,"Symbol Set");
  return null;
}

Tuple!(Individuals, Gref[][]) read_genotype_qtab(File f, ObservedGenotypes symbols) {
  Individuals ret_individuals = new Individuals;
  Gref ret_genotypes[][];  // return matrix
  string buf;
  // note, we skip the marker names - they are for reference only
  while (f.readln(buf)) {
    if (strip(buf) == "# --- Data Observed end")
       break;
    if (buf[0] == '#') continue;
    
    auto res = parse_genotype_qtab(buf);
    auto name_ind_str = res[0];   // string
    auto genotype_ind_str = res[1]; // array of string
    // writeln(individual,"\t",genotype_ind_str);

    // For every genotype symbol (A,B,H,C, etc) we need to find the 
    // matching ObservedGenotype, which is derived from the symbol 
    // table accompanying the genotype data. Next we add the reference
    // to the genotype matrix.
    Gref[] genotype_ind;
    genotype_ind.reserve(genotype_ind_str.length);

    foreach(g ; genotype_ind_str) {
      auto genotype = symbols.decode(g); // genotype is now a reference
      genotype_ind ~= genotype;
    }
    
    ret_individuals.list ~= new Individual(name_ind_str);
    ret_genotypes ~= genotype_ind;
  }
  return tuple(ret_individuals, ret_genotypes);
}

unittest {
  // Symbol and genotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));

  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
  // First read symbol information (the GenotypeCombinators)
  writeln("reading ",symbol_fn);
  auto f = File(symbol_fn,"r");
  // auto symbol_settings = read_set_symbol_qtab(f);
  auto symbols = read_genotype_symbol_qtab(f, false);
  // Test working of symbols
  // assert(symbols.decode("A") == symbols.decode("AA"));
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("H")) == "[(0,1), (1,0)]");
  assert(to!string(symbols.decode("HorA")) == "[(0,0), (0,1), (1,0)]", to!string(symbols.decode("HorA")));
  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"regression","test_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto f1 = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(f1, symbols);
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
}



