/** 
 * This module contains the primitive objects used for QTL mapping.
 *
 * Test the module with 
 *
 *   rdmd --main -unittest qtl/core/primitives.d 
 */

module qtl.core.primitives;

import std.container; 

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
  const uint id;            /// Unique identifier (maybe we don't need this)
  const string name;        /// Name
  Attribute[] attrib_list;  /// Ref. to list of attributes
}

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

struct Marker {
  mixin PayLoad;
  Chromosome chromosome;
  double position;          /// Marker position - content depends on map
}

/**
 * Genotype is the most primitive representation of a genotype. The type
 * can be any type T (normally char or uint, but other objects may be
 * possible).
 *
 * Note the primitive should be small, there may be many genotypes!
 */

struct Genotype(T) {
  T value;
}

/**
 * Phenotype is the most primitive representation of a phenotype. The type
 * can be any type T (normally a double, but can potentially be any Object).
 *
 * Note the primitive should be small, there may be many phenotypes!
 */

struct Phenotype(T) {
  T value;
}

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
  this(string _name, uint _id = -1) {
    id = _id;
    name = _name;
  }
}

class Autosome : Chromosome {
  static is_sex = false;
  this(string _name) { super(_name); };
}

class SexChromosome : Chromosome {
  static is_sex = true;
  this(string _name) { super(_name,0); assert(id == 0); };
}

/**
 * Individual is the most primitive representation of an individual
 * in the context of QTL mapping. Tracking phenotypes etc. is handled
 * by a shared object.
 */

class Individual {
  mixin PayLoad;
}

/******************************************************************************
 * The following objects are not really primitive - but are useful
 * building blocks.
 */

/**
 * The MarkerIndex combines a Marker with a genotype matrix, where each row
 * represents the genotype of an individual.  Each MarkerRef points to a
 * Marker, a genotype matrix, and the column index in the matrix.
 */

struct MarkerRef(T) {
  Marker marker;
  Genotype!T[][] genotype_matrix;
  uint column;
}

/**
 * The ordered Marker list keeps track of MarkerRefs.
 */

class Markers(T) {
  MarkerRef!T[] marker_list;  // Will probably become a List.
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


/******************************************************************************
 * Unit tests for primitives 
 */

import std.stdio;

unittest {
  writeln("Unit test " ~ __FILE__);
  // test marker
  Marker m1 = { id:1, position:4.6};
  assert(m1.id == 1);
  assert(m1.attrib_list == null);
  assert(m1.attrib_list.length == 0);
  Marker m2 = { id:2, position:4.8, chromosome:new Autosome("1"), attrib_list:new Attribute[1]};
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

  // Chromosome
  auto c1 = new Autosome("1");
  auto x = new SexChromosome("X");
  assert(c1.is_sex == false);
  assert(x.is_sex == true);
}

unittest {
  // this should compile 
  // auto map = new FullMap!F2();
  auto map = new FullMap!uint();
  foreach ( c ; map.chromosome_map ) {
    auto markers = c.markers;
    foreach ( m ; markers.marker_list ) {
    }
  }
}
