/**
 * This module contains the primitive objects used for QTL mapping.
 *
 * Test the module with 
 *
 *   rdmd --main -unittest qtl/core/primitives.d 
 */

module qtl.core.primitives;

import std.container; 
import std.conv;
import std.range;
import std.algorithm;
import std.string;
import std.typecons; 

// Default values for undefined types
immutable MARKER_POSITION_UNKNOWN = double.nan;
immutable NAME_UNKNOWN = "unknown";
immutable MARKER_NAME_UNKNOWN = "marker name unknown";
immutable INDIVIDUAL_NAME_UNKNOWN = "individual name unknown";
immutable ID_UNKNOWN = uint.max;

/** 
 * Attribute is a container for additional information that is not
 * anticipated in primitive objects. I.e. Gene Ontology annotation could
 * be tied against a primitive object through the attribute list.
 *
 * Derive your attributes from this class.
 */

class Attribute {
  private string description = "unknown attribute";
}

/** 
 * Primitives have an id, an optional name, and an attribute list.
 */

mixin template Identity()
{
  const uint id;            /// Unique identifier (maybe we don't need this as 
                            /// the memory address is also a unique number)
  const string name;        /// Name
  this() { id = ID_UNKNOWN; name = NAME_UNKNOWN; }
  this(uint _id) { id = _id; name = to!string(_id); }
  this(string _name, uint _id = ID_UNKNOWN) { id = _id; name = _name; }
}

mixin template Attributes()
{
  Attribute[] attrib_list;  /// Ref. to list of attributes
}

// Named primitives - for descriptive parameters
alias double Position;
alias string[string] Founders;

mixin template MarkerInfo() {
  mixin Identity;
  mixin Attributes;
  Chromosome chromosome;      /// Reference to Chromosome
  Position position;          /// Marker position - content depends on map
}

/**
 * Mixing in ActList will create a list of T, and turn the class into an
 * iterator for foreach
 */

mixin template ActList(T) {
  T[] list;
  const size_t length() { return list.length; }

  // The functional way of creating an iterator in D - called by 'foreach'
  int opApply(int delegate(ref T) call) {
    foreach (item ; list) {
      auto result = call(item);
      if (result) return result; // result <> 0 when 'foreach' bails out early
    }
    return 0;
  }

  T add(T item) { 
    list ~= item;
    return item;
  }
}

/** 
 * The Marker struct is the most primitive representation of a marker, i.e.
 * the marker ID, the position, a reference to the marker name, and a reference
 * to a list of attributes - where an attribute can be any object. The
 * attrib_list gives a flexible way of tracking state. Arguably chromosome and 
 * name are attributes, but they are considered to be standard.
 *
 * The Marker object does not keep track of the parent container (apart from
 * Chromosome). Normally a Marker list is maintained in a parent container.
 *
 * Markers are compared on their position (location). I.e. markers at exactly
 * the same location are considered equal (currently the chromosome is not
 * checked).
 */

class Marker {
  mixin MarkerInfo;

  this(double _position = MARKER_POSITION_UNKNOWN, string _name = MARKER_NAME_UNKNOWN, uint _id=ID_UNKNOWN) { 
    name = strip(_name), position = _position, chromosome = null, id = _id ;
    if (name == MARKER_NAME_UNKNOWN && position != MARKER_POSITION_UNKNOWN) 
      // default to 'loc\d\d\d\.\d' name
      name = "loc" ~ to!string(position);
  }
  this(in Marker m) {
    this(m.position, m.name, m.id);
    chromosome = cast(Chromosome)m.chromosome;
  }
  this(Chromosome _chromosome, double _position, string _name = null) {
    chromosome = _chromosome, position = _position;
    id = ID_UNKNOWN;
    name = (_name ? _name : MARKER_NAME_UNKNOWN); // FIXME: auto init
  }
  this(Chromosome _chromosome, double _position, string _name = null, uint _id=ID_UNKNOWN) {
    chromosome = _chromosome, position = _position;
    id = _id;
    name = (_name ? _name : MARKER_NAME_UNKNOWN); // FIXME: auto init
  }

  const Position get_position() { return position; }
  override bool opEquals(Object other) {
    //  test for chromosome?
    /// Markers at the same position are considered equal
    return get_position() == (cast(Marker)other).position;
  }
  override int opCmp(Object other) {
    //  test for chromosome?
    auto cmp = position - (cast(Marker)other).position;
    if (cmp > 0) return 1;
    if (cmp < 0) return -1;
    return 0;
  }
  override string toString() { return name ~ "~" ~ to!string(position); }
}

/** 
 * PseudoMarker is a special marker
 */

class PseudoMarker : Marker {

  this(double _position = MARKER_POSITION_UNKNOWN, string _name = MARKER_NAME_UNKNOWN, uint _id=ID_UNKNOWN) { 
    super(_position, _name, _id);
  }

  this(in PseudoMarker m) {
    this(m.position, m.name, m.id);
  }
  this(Chromosome _chromosome, double _position) 
  { super(_chromosome, _position); }
}

bool isPseudoMarker(M)(in M a)
{
  return(typeid(a) == typeid(PseudoMarker));
}

immutable GENOTYPE_NA = -1; // don't use

/**
 * Note: Genotype is deprecated - use GenotypeCombinator instead!
 * (see genotype.d)
 *
 * Genotype is the most primitive representation of a genotype. The type
 * can be any type T (normally char or uint, but other objects may be
 * possible).
 *
 * Note the primitive should be as small as possible, there may be many 
 * genotypes! Therefore it is a struct (where value may be an int index)
 */

struct Genotype(T) {
  // deprecated
  T value;
  
  /// String representation of genotype.
  string toString() {
    static if (is(typeof(to!int(value)==0)==bool)) {
      // ---- handling when T resolves to an int (or enum). There
      //      is another way - by defining int.GENOTYPE_NA. But
      //      this works, and delegates.
      return to!int(value) == GENOTYPE_NA ? "-" : to!string(value);
    } 
    else {
      // ---- otherwise delegate to genotype implementation
      return to!string(value);
    }
  }
}

alias double Probability;


/**
 * GenoProb keeps track of genotype probabilities
 */

struct GenoProb {
  double value;
  /// String representation of genotype probability.
  string toString(){
    if(to!int(value) != GENOTYPE_NA){
      return to!string(value);
    }else{
      return "-";
    }
  }
}

/**
 * GenoProbs is a three dimensional array holding markers (x-axis),
 * individuals (y-axis) and genotypes (z-axis)
 */

alias GenoProb[][][] GenoProbs;  // = new double[][][](n_gen,n_ind,n_markers);

immutable PHENOTYPE_NA = double.max; // FIXME: needs to be typed to T

/**
 * Phenotype is the most primitive representation of a phenotype. The type
 * can be any type T (normally a double, but can potentially be any Object).
 *
 * Note the primitive should be small as small as possible, there may be many
 * phenotypes! Therefore it is a struct.
 */

struct Phenotype(T) {
  T value;
  
  /// String representation of phenotype.
  const string toString(){
    if(to!double(value) != PHENOTYPE_NA){
      return to!string(value);
    }else{
      return "NA";
    }
  }
}

/**
 * PhenotypeMatrix holds phenotypes (cols) against individuals (rows)
 */

mixin template RealizePhenotypeMatrix(T)
{
  alias Phenotype!T[][] PhenotypeMatrix; // = new double[][][](n_ind,n_phe);
}


/**
 * Covariate representation
 */

struct Covariate(T) {
  T value;

  /// String representation of phenotype.
  string toString(){
    if(value != PHENOTYPE_NA){
      return to!string(value);
    }else{
      return "NA";
    }
  }
}


/**
 * Chromosome is the most primitive representation of a chromosome.  Autosome
 * and sex chromosomes are known via their type. Since these chromosomes are
 * 'shared', or referenced, between markers, we use them by reference 
 * (i.e. a D class).
 *
 * To maintain a list of markers with a chromosome, use a shared object,
 * for example ChromosomeMarkers below.
 */

class Chromosome {
  mixin Identity;

  /// functions, such as is_sex, are defined in chromosome.d
}

class Autosome : Chromosome {
  this(string _name, uint _id) { super(_name,_id); }
}

class SexChromosome : Chromosome {
  this(string _name, uint _id=ID_UNKNOWN) { super(_name,_id); };
}

/**
 * Individual is the most primitive representation of an individual
 * in the context of QTL mapping. Tracking phenotypes etc. is handled
 * by a shared object.
 */

class Individual {
  mixin Identity;

  /// We don't store chromosomes and markers in this object.
  /// See also Markers!M and chromosome.d to 
}

alias string[] Inds;

@property M[] list(M)(M[] ms) { return ms; };

/******************************************************************************
 * Unit tests for primitives 
 */

import std.stdio;

unittest {
  writeln("Unit test " ~ __FILE__);
  // test marker
  Marker m1 = new Marker(4.6,"m1",1);
  assert(m1.id == 1);
  assert(m1.attrib_list == null);
  assert(m1.attrib_list.length == 0);
  Marker m2 = new Marker(4.8);
  m2.chromosome = new Autosome("1",1);
  m2.attrib_list = new Attribute[1];
  assert(m2.attrib_list.length == 1);
  // create a list of markers
  auto markers = [ m1 ];
  markers ~= m2 ;
  assert(markers.length == 2);

  // Genotype
  Genotype!char g1 = { value:'A' };
  assert(g1.value == 'A');
  assert(g1.value != 2);

  // Phenotype
  Phenotype!double p1 = { value:-7.809 };
  assert(p1.value == -7.809);
}

unittest {
  // Test list of markers and pseudomarkers
  Marker m1 = new Marker(4.6,"m1",1);
  Marker m2 = new Marker(4.8,"m2",2);
  m2.chromosome = new Autosome("1",1);
  m2.attrib_list = new Attribute[1];
  PseudoMarker pm1 = new PseudoMarker(4.7,"pm1",3);
}

