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

alias double Position;

mixin template MarkerInfo() {
  mixin Identity;
  mixin Attributes;
  Chromosome chromosome;      /// Reference to Chromosome
  Position position;          /// Marker position - content depends on map
}

mixin template ActList(T) {
  T[] list;
  const uint length() { return list.length; }
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
  this(Chromosome _chromosome, double _position) {
    chromosome = _chromosome, position = _position;
    id = ID_UNKNOWN, name = MARKER_NAME_UNKNOWN; // FIXME: auto init
  }

  Position get_position() { return position; }
  bool opEquals(Object other) {
    //  test for chromosome?
    /// Markers at the same position are considered equal
    return get_position() == (cast(Marker)other).position;
  }
  int opCmp(Object other) {
    //  test for chromosome?
    auto cmp = position - (cast(Marker)other).position;
    if (cmp > 0) return 1;
    if (cmp < 0) return -1;
    return 0;
  }
  string toString() { return name ~ "~" ~ to!string(position); }
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

immutable GENOTYPE_NA = -1;

/**
 * Genotype is the most primitive representation of a genotype. The type
 * can be any type T (normally char or uint, but other objects may be
 * possible).
 *
 * Note the primitive should be as small as possible, there may be many 
 * genotypes! Therefore it is a struct.
 */

struct Genotype(T) {
  T value;
  
  /// String representation of genotype.
  string toString(){
    if(to!int(value) != GENOTYPE_NA){
      return to!string(value);
    }else{
      return "-";
    }
  }
}

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
  string toString(){
    if(to!double(value) != PHENOTYPE_NA){
      return to!string(value);
    }else{
      return "NA";
    }
  }
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

class Individuals {
  mixin ActList!Individual;
}

/******************************************************************************
 * The following objects are not really primitive - but are the
 * building blocks tying primitive types together.
 */

/**
 * Combine a Marker with a genotype matrix. Each MarkerRef points to a Marker,
 * a genotype matrix, and the column index in that matrix. 
 *
 * Different MarkerRefs can point to different genotype matrices(!) E.g. a
 * special genotype matrix for one or more pseudomarkers can exist.
 */

class MarkerRef(T) 
{
  Marker marker;
  Genotype!T[][] genotype_matrix;
  uint column;

  this(double _position, string _name = MARKER_NAME_UNKNOWN, uint _id=ID_UNKNOWN) { 
    marker = new Marker(_position, _name, _id); }
  this(in MarkerRef!T mr) {
    this(mr.marker);
  }
  this(in Marker m) {
    marker = new Marker(m.position, m.name, m.id);
  }

  Position get_position() { return marker.get_position; }
  const string name() { return marker.name; }
  const uint id() { return marker.id; }

  /// Markers at the same position are considered equal
  bool opEquals(Object other) {
    return get_position() == (cast(MarkerRef!T)other).get_position();
  }
  int opCmp(Object other) {
    auto cmp = get_position() - (cast(MarkerRef!T)other).get_position();
    if (cmp > 0) return 1;
    if (cmp < 0) return -1;
    return 0;
  }

}

/**
 * The Marker list keeps track of Markers. Note there is no 
 * guarantee the list is ordered.
 */

class Markers(M) {
  M[] list;  // Unordered marker list May become an SList.
  this() {}
  this(in Markers!M markers) {
    list = markers.list.dup;  // make sure to clone all data
  }
  this(M[] markers) {
    list = markers.dup;
  }
  void add(in Marker m) {
    list ~= new M(m);
  }
  const auto sorted() { 
    auto ms = new Markers(this); // make a copy
    sort(ms.list); // sorts in place
    return ms;
  }
  const int length() { return list.length; }
  M opIndex(int i) { return list[i]; }
  /// find by name
  M find(string name) {
    foreach(m; list) { if (m.name == name) return m; }
    return null;
  }
  /// find by ID
  M find(uint id) {
    foreach(m; list) { if (m.id == id) return m; }
    return null;
  }
}

@property M[] list(M)(Markers!M ms) { return ms.list; }
@property M[] list(M)(M[] ms) { return ms; };


/**
 * ChromosomeMap combines Chromosome and Marker list.
 */

class ChromosomeMap(T) {
  Chromosome chromosome;
  // Markers!(MarkerRef!T) markers;
}

/**
 * The FullMap has an ordered chromosome map.
 */

class FullMap(T) {
  SList!(ChromosomeMap!T) chromosome_map;
}

class MappedQTLs {
}

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
  // e.g. for F2 auto map = new FullMap!F2();
  // this should also compile:
  auto map = new FullMap!uint();
  foreach ( c ; map.chromosome_map ) {
  /*
    auto markers = c.markers;
    foreach ( m ; markers.list ) {
    }
  */
  }
}

unittest {
  // Test list of markers and pseudomarkers
  Marker m1 = new Marker(4.6,"m1",1);
  Marker m2 = new Marker(4.8,"m2",2);
  m2.chromosome = new Autosome("1",1);
  m2.attrib_list = new Attribute[1];
  auto mref1 = new MarkerRef!uint(m1);
  auto mref2 = new MarkerRef!uint(m2);
  PseudoMarker pm1 = new PseudoMarker(4.7,"pm1",3);
  auto pmref1 = new MarkerRef!uint(pm1);
  auto markers = new Markers!(MarkerRef!uint)();
  markers.list ~= mref1;  // 1, 4.6
  markers.list ~= mref2;  // 2, 4.8
  markers.list ~= pmref1; // 3, 4.7
  // test list
  assert(list(markers) == markers.list);
  // find by index
  assert(markers[0].name == "m1");
  // find by name
  assert(markers.find("m1").name == "m1");
  uint[] result;
  foreach ( m ; markers.list ) {
    result ~= m.marker.id;
  }
  assert(result==cast(uint[])[1,2,3]);
  // reverse sort markers
  auto list2 = sort!("a.get_position() > b.get_position()")(markers.list);
  foreach ( m ; list2 ) {
    result ~= m.marker.id; // add reverse sorted list (2,3,1)

  }
  assert(result==cast(uint[])[1,2,3,2,3,1],to!string(result));
  auto ms = markers.sorted(); 
  foreach ( m ; ms.list ) {
    result ~= m.marker.id; // add sorted marker list (1,3,2)
  }
  assert(result==cast(uint[])[1,2,3,2,3,1,1,3,2],to!string(result));
  auto ulist2 = uniq!("a.get_position() == b.get_position")(ms.list);
  auto pos_list = map!"a.get_position()"(ulist2);
  assert(equal(pos_list,[4.6, 4.7, 4.8]),to!string(pos_list));
}
