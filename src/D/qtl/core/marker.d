/**
 * Marker related functions
 */


module qtl.core.marker

import std.container;
import std.conv;
import std.variant;
import std.random;
import qtl.core.primitives;
import std.algorithm; // to use sort()

import std.stdio;
// import std.container;
import std.typecons;

/**
 * The Marker list keeps track of Markers. Note there is no 
 * guarantee the list is ordered.
 */

class Markers(M) {
  mixin ActList!M; 
  this() {}
  this(in Markers!M markers) {
    // list = markers.list.dup;  // make sure to clone all data
    foreach(m ; markers.list) {
      list ~= cast(M) m;
    }
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
  const size_t length() { return list.length; }
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

Marker[] dup(in Marker[] ms) 
{ 
  Marker[] new_ms;
  foreach(m; ms) {
    new_ms ~= cast(Marker)m;
  }
  return new_ms;
}



