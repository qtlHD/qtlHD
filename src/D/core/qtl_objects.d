/** This module contains the primitive objects used for QTL mapping.
 */

module qtl.object;

/** 
 * Attribute is a container for additional information that is not
 * anticipated in primitive objects. I.e. Gene Ontology annotation could
 * be tied against a primitive object through the attribute list.
 *
 * Derive your attributes from this class.
 */

class Attribute {
}

/** 
 * The Marker struct is the most primitive representation of a marker, i.e.
 * the marker ID, the position, a reference to the marker name, and a reference
 * to a list of attributes - where an attribute can be any object. The
 * attrib_list gives a flexible way of tracking state.
 *
 * The Marker object does not auto
 */

struct Marker {
  const uint id;            /// Unique identifier
  const int chromosome;
  double position;   /// The marker position - content depends on map
  Attribute[] attrib_list;
  unittest {
    Marker m = { id:1, position=4.6};	      
  }
}


