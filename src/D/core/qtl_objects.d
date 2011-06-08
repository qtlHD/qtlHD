/** 
 * This module contains the primitive objects used for QTL mapping.
 *
 * Test the module with 'rdmd --main -unittest qtl_objects.d'
 */

module qtl.qtl_object;

/** 
 * Attribute is a container for additional information that is not
 * anticipated in primitive objects. I.e. Gene Ontology annotation could
 * be tied against a primitive object through the attribute list.
 *
 * Derive your attributes from this class.
 */

class Attribute {
  private string description = "core attribute";
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
  const uint id;            /// Unique identifier
  const string name;        /// Marker name
  const int chromosome;
  double position;          /// Marker position - content depends on map
  Attribute[] attrib_list;  /// Ref. to list of attributes
}

/**
 * Chromosome is the most primitive representation of a chromosome. 
 * Sex chromosomes are known via their name.
 */

struct Chromosome {
  const uint id;            /// Unique identifier
  const string name;        /// Chromosome name
  Attribute[] attrib_list;  /// Ref. to list of attributes
}

/**
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
  // create a list
  auto markers = [ m1 ];
  markers ~= m2 ;
  assert(markers.length == 2);
}

