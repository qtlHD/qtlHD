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
import std.math; 

// Default values for undefined types
immutable VALUE_NAN = double.nan;
immutable MARKER_POSITION_UNKNOWN = VALUE_NAN;
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
  immutable uint id;        /// Unique identifier (maybe we don't need this as 
                            /// the memory address is also a unique number)
  immutable string name;    /// Name
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

  /*
  Individuals opOpAssign(string op)(string name) if (op == "~") {
    list ~= new Individual(name);
    return this;
  }

  Individuals opOpAssign(string op)(Individual ind) if (op == "~") {
    // list ~= ind.dup;
    return this;
  }
  */
 }

/**
 * Class for supporting Range iterators on classes that contain
 * a list (e.g. for using map!)
 */

class ListIter(T) {

  uint cursor = 0;
  T contained;

  this (T container) {
    contained = container;
  }

  @property auto front() {
    // writeln("Access ",cursor);
    return contained.list[cursor];
  }
 
  void popFront() {
    // writeln("Cursor ",cursor);
    cursor += 1;
  }
 
  @property bool empty() { 
    // writeln("Cursor,len ",cursor,",",contained.list.length);
    return cursor >= contained.list.length; 
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

alias double Probability;

/**
 * Covariate representation
 */

struct Covariate(T) {
  T value;

  string toString(){
    if(isNaN(value)){
      return to!string(value);
    }else{
      return "NA";
    }
  }
}


/**
 * Chromosome is the most primitive representation of a chromosome.  Autosome
 * and X chromosomes are known via their type. Since these chromosomes are
 * 'shared', or referenced, between markers, we use them by reference 
 * (i.e. a D class).
 *
 * To maintain a list of markers with a chromosome, use a shared object,
 * for example ChromosomeMarkers below.
 */

class Chromosome {
  mixin Identity;

  /// functions, such as is_X_chr, are defined in chromosome.d
}

class Autosome : Chromosome {
  this(string _name, uint _id) { super(_name,_id); }
}

class Xchromosome : Chromosome {
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


