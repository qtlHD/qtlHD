/** 
 * This module contains the primitive objects used for QTL mapping.
 *
 * Test the module with 
 *
 *   rdmd --main -unittest qtl/core/primitives.d 
 */

module qtl.core.primitives;

import std.container; 

immutable MARKER_POSITION_UNKNOWN = double.nan;
immutable MARKER_NAME_UNKNOWN = "unknown";
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

mixin template PayLoad()
{
  const uint id;            /// Unique identifier (maybe we don't need this as 
                            /// the memory address is also a unique number)
  const string name;        /// Name
  Attribute[] attrib_list;  /// Ref. to list of attributes
}

alias double Position;

/** 
 * The Marker struct is the most primitive representation of a marker, i.e.
 * the marker ID, the position, a reference to the marker name, and a reference
 * to a list of attributes - where an attribute can be any object. The
 * attrib_list gives a flexible way of tracking state. Arguably chromosome and 
 * name are attributes, but they are pretty standard.
 *
 * The Marker object does not keep track of the parent container. Normally
 * a Marker list is maintained in a parent container.
 */

class Marker {
  mixin PayLoad;
  bool is_pseudo() { return false; }; // I would like to get this from the type system
  Chromosome chromosome;
  Position position;          /// Marker position - content depends on map
  this(double _position = MARKER_POSITION_UNKNOWN, uint _id=ID_UNKNOWN, string _name = MARKER_NAME_UNKNOWN) { 
    name = _name, position = _position, id = _id ;
  }
}

class PseudoMarker : Marker {
  override bool is_pseudo() { return true; };
  this(double _position = MARKER_POSITION_UNKNOWN, uint _id=0, string _name = "unknown") { super(_position, _id, _name); }
}

/**
 * Genotype is the most primitive representation of a genotype. The type
 * can be any type T (normally char or uint, but other objects may be
 * possible).
 *
 * Note the primitive should be small, there may be many genotypes! Therefore
 * it is a struct.
 */

struct Genotype(T) {
  T value;
}

/**
 * Phenotype is the most primitive representation of a phenotype. The type
 * can be any type T (normally a double, but can potentially be any Object).
 *
 * Note the primitive should be small, there may be many phenotypes! Therefore 
 * it is a struct.
 */

struct Phenotype(T) {
  T value;
}

/** 
 * Phenotype objects.
 *
 *

interface PhenotypeContainer {
}

class Phenotypes(T) : PhenotypeContainer {
  Phenotype!T[] list;
}
*/

struct Covariate(T) {
  T value;
}

/*
interface CovariateContainer {
}

class Covariates(T) : CovariateContainer {
  Covariate!T[] list;
}

*/

import std.conv;

/**
 * Chromosome is the most primitive representation of a chromosome.
 * Autosome and sex chromosomes are known via their type. Since these
 * chromoses are 'shared' between markers, we use them by reference (i.e. a
 * class).  
 *
 * To maintain a list of markers with a chromosome, use a shared object,
 * for example ChromosomeMarkers below.
 */

class Chromosome {
  mixin PayLoad;
  bool is_sex() { return false; };
  this(string _name, uint _id = -1) {
    id = _id;
    name = _name;
  }
}

class Autosome : Chromosome {
  this(string _name, uint _id) { super(_name,_id); assert(!is_sex); }
}

class SexChromosome : Chromosome {
  override bool is_sex() { return true; };
  this(string _name, uint _id=0) { super(_name,_id); assert(is_sex); };
}

/**
 * Individual is the most primitive representation of an individual
 * in the context of QTL mapping. Tracking phenotypes etc. is handled
 * by a shared object.
 */

class Individual {
  mixin PayLoad;
}

class Individuals {
  SList!Individual list;
}

/******************************************************************************
 * The following objects are not really primitive - but are useful
 * building blocks.
 */

/**
 * The MarkerIndex combines a Marker with a genotype matrix, where each row
 * represents the genotype of an individual.  Each MarkerRef points to a
 * Marker, a genotype matrix, and the column index in the matrix. 
 *
 * Multiple genotype matrices are possible. E.g. a genotype matrix for one or 
 * more pseudomarkers can exist.
 */

class MarkerRef(T) {
  Marker marker;
  Genotype!T[][] genotype_matrix;
  uint column;
  this(double _position, uint _id=ID_UNKNOWN, string _name = MARKER_NAME_UNKNOWN) { 
    marker = new Marker(_position, _id, _name); }
  this(double _position, string _name) { 
    marker = new Marker(_position, ID_UNKNOWN, _name);
  }
  this(Marker m) {
    marker = m;
  }
}

/**
 * The ordered Marker list keeps track of MarkerRefs.
 */

interface MarkerContainer {
}

class Markers(T) : MarkerContainer {
  MarkerRef!T[] list;  // Will probably become a List.
  auto markercontainer() { return list; }
  this() {}
  this(in MarkerRef!T[] _list) {
    list = _list.dup;
  }
}

/**
 * ChromosomeMap combines Chromosome and Marker list.
 */

class ChromosomeMap(T) {
  Chromosome chromosome;
  Markers!T markers;
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
  Marker m1 = new Marker(4.6,1);
  assert(m1.id == 1);
  assert(m1.attrib_list == null);
  assert(m1.attrib_list.length == 0);
  Marker m2 = new Marker(4.8,2);
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
    auto markers = c.markers;
    foreach ( m ; markers.list ) {
    }
  }
}

unittest {
  // Test list of markers and pseudomarkers
  Marker m1 = new Marker(4.6,1);
  Marker m2 = new Marker(4.8,2);
  m2.chromosome = new Autosome("1",1);
  m2.attrib_list = new Attribute[1];
  auto mref1 = new MarkerRef!uint(m1);
  auto mref2 = new MarkerRef!uint(m2);
  PseudoMarker pm1 = new PseudoMarker(4.7,3);
  auto pmref1 = new MarkerRef!uint(pm1);
  auto markers = new Markers!uint();
  markers.list ~= mref1;
  markers.list ~= mref2;
  markers.list ~= pmref1;
  uint[] result;
  foreach ( m ; markers.list ) {
    result ~= m.marker.id;
  }
  assert(result==cast(uint[])[1,2,3]);
}
