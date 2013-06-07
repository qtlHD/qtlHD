/**
 * Read CSVR format
 *
 * Note: This implementation handles everything through one object. This may change later.
 **/

module qtl.plugins.csvr.read_csvr;

import core.memory;

import core.memory;

import qtl.core.primitives;
import qtl.core.chromosome;
// import qtl.core.util.matrix;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.individual;
import example.genotype_examples;

import std.stdio;
import std.conv;
import std.string;
import std.path;

/** 
 * Loads a simple CSVR file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the multitrait.csvr file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */
class CSVrReader(XType,ObservedXType) {
  private File f;
  XType crosstype;
  ObservedXType symbols;
  string[] phenotypenames;
  Marker[] markers;
  Individuals individuals;
  Chromosome[string] chromosomes;
  Phenotype[][] phenotypes;
  GenotypeCombinator[][] genotypes;
  size_t n_phenotypes;
  
  bool is_phenotype(string[] location){ return(location == ["",""]); }
  
  this(in string fn, ObservedXType observed = null){
    f = File(fn,"r");
    scope(exit) f.close();
    
    if (observed is null){
      symbols = new ObservedXType;
    }else{
      symbols = observed;
    }
    writeln("SYMBOLS:",symbols);
    
    string line;
    int linecount;
    while((line = f.readln()) != ""){
      line = strip(line);
      debug writeln("LINE:",line);
      auto fields = split(line,",");

      if(is_phenotype(fields[1..3])){
        //Phenotype
        debug writeln("Phenotype: " ~ fields[0]);
        Phenotype[] ps;
        phenotypenames ~= fields[0];
        foreach (field; fields[3..$]) {
          ps ~= set_phenotype(strip(field));
        }
        phenotypes ~= ps;
      }else{
        debug writeln("Marker: " ~ fields[0]);
        //CHR
        if (!(fields[1] in chromosomes)){
          chromosomes[fields[1]] = get_chromosome_with_id(fields[1]);
        }

        //Marker
        Marker m = new Marker(to!double(fields[2]),fields[0]); //
        m.chromosome = chromosomes[fields[1]];
        markers ~= m;

        //Genotype
        GenotypeCombinator[] gs;
        // we use the predefined crosstype symbols
        foreach (field; fields[3..$]) {
          gs ~= symbols.decode(strip(field));
        }
        genotypes ~= gs;
      }
      linecount++;
    }
    writefln("Read %s with %d phenotypes and %d markers measured at %d individuals",fn, phenotypes.length, genotypes.length, genotypes[0].length);
    f.close();
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.buildPath buildPath;
  auto infn = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..","..","test","data","input","multitrait.csvr"));
  writeln("  - reading CSVR " ~ infn);
  auto data = new CSVrReader!(RIL,ObservedRIL)(infn);
  assert(data.phenotypes.length == 24, to!string(data.phenotypes.length));
  assert(data.markers.length == 117, to!string(data.markers.length));
  assert(data.phenotypes[0].length == 162, to!string(data.phenotypes[0].length));
  assert(data.genotypes[0].length == 162, to!string(data.genotypes[0].length));
}

