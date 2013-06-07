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
import std.variant;
alias std.string.split str_split;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.individual;

// mixin RealizePhenotypeMatrix!double;

// alias Phenotype[][] PhenotypeMatrix;

/**
 * Low level parsers, str_splits a tab delimited line into fields,
 * and removes everything after a remark '#' symbol. Optionally
 * strips each field of white space.
 */


string[] split_line(string line,string split_with="\t") {
  // strip remark and leading whitespace
  auto res = splitter(chomp(line), regex("\\s+#\\s"));
  return str_split(to!string(std.array.array(res)[0]),split_with);
}

string[] split_line_on_spaces(string line) {
  return split_line(line," ");
}

string[] split_line_strip(string line) {
  auto fields1 = split_line(line);
  return std.array.array(map!"strip(a)"(fields1));
}

string[] split_line_on_whitespace(string line) {
  auto res = splitter(strip(line), regex("\\s+#\\s"));
  auto content = to!string(std.array.array(res)[0]);
  auto res1 = splitter(content, regex("\\s+"));
  return std.array.array(res1);
}

unittest {
  writeln("Unit test " ~ __FILE__);
  assert(split_line("test  # remark") == ["test"],to!string(split_line("test  # remark")));
  assert(split_line("test\n") == ["test"],to!string(split_line("test")));
  assert(split_line("test\ttest") == ["test","test"]);
  assert(split_line("test\t test \t ") == ["test"," test "," "]);
  assert(split_line_strip("test\t test\t") == ["test","test",""]);
  assert(split_line_strip("test\t test\t") == ["test","test",""]);
  assert(split_line_strip("test\t test\t# remark") == ["test","test"]);
  assert(split_line_strip("test\t test\t # remark") == ["test","test"]);
  assert(split_line_on_spaces("test test\t # remark") == ["test","test"]);
  writeln(split_line_on_spaces("test test  test \n"));
  assert(split_line_on_spaces("test test  test \n") == ["test","test","","test",""]);
  writeln(split_line_on_whitespace("test test\t test  test \n"));
  assert(split_line_on_whitespace("test test\t test  test \n") == ["test","test","test","test"]);
}

/** 
 * Parse a line for key values
 */

Tuple!(string, string[]) parse_line_key_values(string line) {
  auto fields = split_line_strip(line);
  auto key = (fields.length > 0 ? fields[0] : null);
  auto value = (fields.length > 1 ? fields[1..$] : null);
  return tuple(key,value);
}

Tuple!(string, string[]) parse_line_key_values_on_whitespace(string line) {
  auto fields = split_line_on_whitespace(line);
  auto key = (fields.length > 0 ? fields[0] : null);
  auto value = (fields.length > 1 ? fields[1..$] : null);
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
        call_line(buf);
      }
    }
  }
}

void each_line_in_section(string fn, string tag, void delegate (string) call_line) {
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  each_line_in_section(f,tag,call_line);
}

/**
 * Turn a QTAB file section into key-value pairs
 */

string[string] get_section_key_values(File f, string tag) {
  string[string] ret;
  each_line_in_section(f,tag, 
  (string line) {
    auto res = parse_line_key_values_on_whitespace(line);
    // writeln(res);
    ret[res[0]] = res[1][0]; // note: only one value per key, this may change
  });
  return ret;
}

string[string] get_section_key_values(string fn, string tag) {
  string[string] ret;
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  return get_section_key_values(f,tag);
}

/**
 * Iterate for section key values without loading the full table in memory twice
 */
void each_section_key_values(string fn, string tag, void delegate (string, string[]) call_line) {
  each_line_in_section(fn, tag, 
    (line) {
      auto res = parse_line_key_values_on_whitespace(line);
      call_line(res[0],res[1]);
    }
  );
}

string[string] read_founder_settings_qtab(string founder_fn) {
  return get_section_key_values(founder_fn,"Set Founder");
}

unittest {
  // Founder reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));

  auto founder_fn = to!string(buildPath(dir,"input","listeria_qtab","listeria_founder.qtab"));
  writeln("reading ",founder_fn);
  auto info = get_section_key_values(founder_fn,"Set Founder");
  writeln(info);
  assert(info["Cross"] == "F2");
}

Tuple!(string, string, double) parse_marker_qtab(string line) {
  auto res = parse_line_key_values(line);
  auto name = res[0];
  auto fields = res[1];
  auto chromosome = (fields.length > 0 ? fields[0] : null);
  auto position = (fields.length > 1 ? to!double(fields[1]) : MARKER_POSITION_UNKNOWN);
  return tuple(name,chromosome,position);
}


Ms[] read_marker_map_qtab(Ms)(string fn) {  // Ms is Marker[] (vs Markers)
  Ms ret_ms[];
  Chromosome[string] clist; // track chromosome objects
  uint id=0;
  each_line_in_section(fn,"Data Location",
    (string line) {
      auto res = parse_marker_qtab(line);
      auto name = res[0];
      auto cname = res[1];
      auto c = (cname in clist ? clist[cname] : new Chromosome(cname));
      auto pos = res[2];
      auto marker = new Marker(c,pos,name,id);
    	id++;
	    ret_ms ~= marker;
    }
  );
  return ret_ms;
}

/**
 * Low level qtlHD qtab parser. Read the file in sections, and
 * delegate to the appropriate readers. Values are allocated
 * in memory and returned in a Tuple.
 */
Tuple!(P[][]) read_phenotype_qtab(P)(string fn) {
  P ret_phenotypes[][];  // return matrix
  each_line_in_section(fn,"Data Phenotype", 
    (string line) {
      auto res = parse_line_key_values(line);
      auto ind = res[0];
      auto fields  = res[1];
      P[] ps = std.array.array(map!((a) {return set_phenotype(a);})(fields));
      ret_phenotypes ~= ps;
    }
  );
  return tuple(ret_phenotypes);
}

unittest {
  // Phenotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));
  auto pheno_fn = to!string(buildPath(dir,"regression","test_phenotype.qtab"));
  writeln("reading ",pheno_fn);
  auto p_res = read_phenotype_qtab!(Phenotype)(pheno_fn);
  Phenotype[][] pheno = p_res[0];
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
  auto fields = split_line_on_spaces(line);
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
    (string line) {
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

ObservedGenotypes read_genotype_symbol_qtab(string filename, bool phase_known = true) {
  auto f = File(filename, "r");
  return read_genotype_symbol_qtab(f, phase_known);
}

class SymbolSettings {
  bool phase_known = false; // default
  override const string toString() {
    return "Phase known: "~(phase_known?"True":"False");
  }
}

/**
 * Read set symbol section in qtab file and convert to a SetSymbol
 * object
 */

SymbolSettings read_set_symbol_qtab(File f) {
  writeln("read_set_symbol_qtab");
  auto settings = new SymbolSettings;
  auto kvs = get_section_key_values(f,"Set Symbol");
  foreach(k, v ; kvs) {
    writeln("SymbolSetting: ",k,v);
  }
  if ("Phase" in kvs) settings.phase_known = (kvs["Phase"] == "known");
  return settings;
}

Tuple!(Individuals, Gref[][]) read_genotype_qtab(File f, ObservedGenotypes symbols) {
  Individuals ret_individuals = new Individuals;
  Gref ret_genotypes[][];  // return matrix
  each_line_in_section(f,"Data Genotype", 
    (string line) {
      // note, we skip the marker names - they are for reference only
      auto res = parse_line_key_values(line);
      auto name_ind_str = res[0];   // string
      auto genotype_ind_str = res[1]; // array of string
      // writeln("Line: ",name_ind_str,"\t",genotype_ind_str);

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
  );
  return tuple(ret_individuals, ret_genotypes);
}

Tuple!(Individuals, Gref[][]) read_genotype_qtab(string filename, ObservedGenotypes symbols) {
  File f = File(filename, "r");
  return read_genotype_qtab(f, symbols);
}

unittest {
  // Symbol and genotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));

  auto symbol_fn1 = to!string(buildPath(dir,"regression","test_symbol_phase.qtab"));
  writeln("reading ",symbol_fn1);
  auto f1 = File(symbol_fn1,"r");
  auto symbol_settings1 = read_set_symbol_qtab(f1);
  assert(symbol_settings1.phase_known == true);
  auto symbols1 = read_genotype_symbol_qtab(f1,symbol_settings1.phase_known);
  assert(to!string(symbols1.decode("AB")) == "[(0,1)]");
  assert(to!string(symbols1.decode("BA")) == "[(1,0)]");

  // First read symbol information (the GenotypeCombinators)
  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
  writeln("reading ",symbol_fn);
  auto f = File(symbol_fn,"r");
  auto symbol_settings = read_set_symbol_qtab(f);
  assert(symbol_settings.phase_known == false);
  auto symbols = read_genotype_symbol_qtab(f,symbol_settings.phase_known);
  // Test working of symbols
  // assert(symbols.decode("A") == symbols.decode("AA"));
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("H")) == "[(0,1), (1,0)]");
  assert(to!string(symbols.decode("HorA")) == "[(0,0), (0,1), (1,0)]", to!string(symbols.decode("HorA")));
  // Read genotype matrix
  auto genotype_fn = to!string(buildPath(dir,"regression","test_genotype.qtab"));
  writeln("reading ",genotype_fn);
  auto f2 = File(genotype_fn,"r");
  auto ret = read_genotype_qtab(f2, symbols);
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

/**
 * Get phenotype matrix by iterating through the data section and building up
 * the matrix, using the lower level parser. At this point only double is
 * supported. Multiple data types will be stored as multiple tables.
 */

Tuple!(string[],PhenotypeMatrix) get_phenotype_matrix(string fn) {
  PhenotypeMatrix p;
  string[] phenotypenames;
  auto type = get_section_key_values(fn,"Type Phenotype");
  uint i = 0;
  each_section_key_values(fn,"Data Phenotype", 
    (key, values) {
      phenotypenames ~= key;
      Phenotype[] ps;
      ps.reserve(values.length);
      foreach (j, v ; values) {
        ps ~= set_phenotype(v);
      }
      i++;
      p ~= ps;
    }
  );
  return tuple(phenotypenames,p);
}

/**
 * Return genotypes as the original strings
 */

Tuple!(string[],string[][]) get_genotype_matrix(string fn) {
  string[][] gm;
  string[] inds;
  // auto type = get_section_key_values(fn,"Type Genotype");
  each_section_key_values(fn,"Data Genotype", 
    (key, values) {
      inds ~= key;
      gm ~= values;
    }
  );
  return tuple(inds,gm);
}

/**
 * Autodetect the contents of a qtab file - based on the first line descriptor
 */

enum QtabFileType {
  symbols,
  founder,
  genotype,
  phenotype,
  location,  // i.e. marker data
  undefined
};

QtabFileType autodetect_qtab_file_type_from_header(string fn, string line) {
  auto fields = split_line_on_whitespace(line);
  writeln(line);
  if (fields[0] != "#" || fields[1] != "---" || fields.length < 5)
    throw new Exception("Malformed detection line in qtab file "~fn~": "~line);
  writeln(fields);
  switch (fields[3]) {
    case "Symbol":    return QtabFileType.symbols;
    case "Phenotype": return QtabFileType.phenotype;
    case "Genotype":  return QtabFileType.genotype;
    case "Founder":   return QtabFileType.founder;
    case "Location":  return QtabFileType.location;
    default:          
      throw new Exception("Cannot autodetect type from qtab file "~fn~": "~line);
  }
}

QtabFileType autodetect_qtab_file_type(in string fn) {
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  string line;
  if (!f.readln(line))
    throw new Exception("Can not autdetect qtab file "~fn);
  return autodetect_qtab_file_type_from_header(fn,line);
}

unittest {
  writeln("autodetect qtab files");
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));

  assert(autodetect_qtab_file_type_from_header("no file","# --- qtl* Symbol Test")==QtabFileType.symbols);
  assert(autodetect_qtab_file_type_from_header("no file","# --- qtl* Genotype Test")==QtabFileType.genotype);

  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol_phase.qtab"));
  // Auto detect read files
  assert(autodetect_qtab_file_type(symbol_fn)==QtabFileType.symbols);
}

/**
 * Note the use of Variant is a health hazard
 */

Variant[] load_qtab(string fn) {
  auto t = autodetect_qtab_file_type(fn);
  with (QtabFileType) {
    switch (t) {
      case symbols: 
        auto f = File(fn,"r");
        auto symbol_settings = read_set_symbol_qtab(f);
        auto observed_genotypes = read_genotype_symbol_qtab(f,symbol_settings.phase_known);
        return variantArray(t,symbol_settings,observed_genotypes);
      case founder: 
        return variantArray(t,get_section_key_values(fn,"Set Founder"));
      case genotype: 
        auto res = get_genotype_matrix(fn);
        return variantArray(t,res[0],res[1]);
      case phenotype: 
        auto res = get_phenotype_matrix(fn);
        return variantArray(t,res[0],res[1]);
      case location: 
        auto ms = read_marker_map_qtab!Marker(fn);
        return variantArray(t,ms);
      default: return null; // tuple(QtabFileType.undefined,variantArray(null));
    }
  }
}

unittest {
  writeln("automatic data loading");
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data"));

  auto symbols = load_qtab(to!string(buildPath(dir,"regression","test_symbol_phase.qtab")));
  assert(symbols[0]==QtabFileType.symbols);
  auto symbol_settings = symbols[1].get!SymbolSettings;  // Class needs the 'cast'
  assert(symbol_settings.phase_known == true);
  auto founders = load_qtab(to!string(buildPath(dir,"input","listeria_qtab","listeria_founder.qtab")));
  assert(founders[0]==QtabFileType.founder);
  auto genotypes = load_qtab(to!string(buildPath(dir,"input","listeria_qtab","listeria_genotype.qtab")));
  assert(genotypes[0]==QtabFileType.genotype);
  auto phenotypes = load_qtab(to!string(buildPath(dir,"input","listeria_qtab","listeria_phenotype.qtab")));
  assert(phenotypes[0]==QtabFileType.phenotype);
  auto pnames = phenotypes[1];
  assert(to!string(pnames[0]) == "1",to!string(pnames[0]));
  auto pmatrix = phenotypes[2];
  assert(to!string(pmatrix[0][0]) == "118.317",to!string(pmatrix[0][0]));
  auto markermap = load_qtab(to!string(buildPath(dir,"input","listeria_qtab","listeria_marker_map.qtab")));
  assert(markermap[0]==QtabFileType.location);
  // auto markers = markermap[1];
  // assert(markers["D15M34"]=="15");
}

/**
 * Main qtab file parser - should load all sections and return the data in the proper
 * containers.
 */

Tuple!(SymbolSettings, Founders, Marker[], Inds, PhenotypeMatrix, ObservedGenotypes, GenotypeMatrix) load_qtab(string[] fns) {
  SymbolSettings s;
  Founders f;
  Marker[] ms;
  Inds i;
  PhenotypeMatrix p;
  string[][] g;
  ObservedGenotypes observed;
  foreach (fn ; fns) {
    auto res = load_qtab(fn);
    auto t = res[0].get!QtabFileType;
    auto d = res[1];
    with (QtabFileType) {
      switch(t) {
        case symbols: 
          s = d.get!SymbolSettings; 
          observed = res[2].get!ObservedGenotypes;
          break;
        case founder: f = d.get!Founders; writeln(f); break;
        case location: 
          ms = d.get!(Marker[]); 
          writeln(ms);
          break;
        case genotype: 
          i = d.get!Inds;
          g = res[2].get!(string[][]);
          break;
        case phenotype: 
          // auto pids = d; ignored, for now
          p = res[2].get!PhenotypeMatrix; break;
        default: throw new Exception("Unsupported file type for " ~ fn);
      }
    }
  }
  // Turn the genotype matrix into a genotype combinator matrix
  auto gc = convert_to_combinator_matrix(g,observed);
  return tuple(s,f,ms,i,p,observed,gc);
}

unittest {
  writeln("automatic data loading (2)");
}
