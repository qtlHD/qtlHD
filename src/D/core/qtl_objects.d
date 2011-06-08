/** 
 * This module contains the primitive objects used for QTL mapping.
 *
 * Test the module with 
 *
 *   rdmd --main -unittest qtl_objects.d
 */

module qtl.objects;

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

mixin template PrimitiveInfo()
{
  const uint id;            /// Unique identifier
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
  mixin PrimitiveInfo;
  const int chromosome;
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

/**
 * Chromosome is the most primitive representation of a chromosome. 
 * Sex chromosomes are known via their name. 
 *
 * To maintain a list of markers with a chromosome, use a shared object.
 */

struct Chromosome {
  mixin PrimitiveInfo;
}

/**
 * Individual is the most primitive representation of an individual
 * in the context of QTL mapping. Tracking phenotypes etc. is handled
 * by a shared object.
 */

struct Individual {
  mixin PrimitiveInfo;
}

/******************************************************************************
 * Unit tests for primitives 
 */

unittest {
  // test marker
  Marker m1 = { id:1, position:4.6};
  assert(m1.id == 1);
  assert(m1.attrib_list == null);
  assert(m1.attrib_list.length == 0);
  Marker m2 = { id:2, position:4.8, chromosome:1, attrib_list:new Attribute[1]};
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

