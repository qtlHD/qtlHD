/**
 * Individual related functions
 */


module qtl.core.individual;

import std.container;
import std.conv;
import std.variant;
import std.random;
import qtl.core.primitives;
import std.algorithm; // to use sort()

import std.stdio;
// import std.container;
import std.typecons;

// FIXME: move cursor out of class

class Individuals {
  uint cursor = 0;
  Individual[] list;
  // mixin ActList!Individual;
  Individuals opOpAssign(string op)(string name) if (op == "~") {
    list ~= new Individual(name);
    return this;
  }
  Individuals opOpAssign(string op)(Individual ind) if (op == "~") {
    writeln(ind);
    // list ~= ind.dup;
    return this;
  }
 
  @property auto front() {
    writeln("Cursor",cursor);
    return list[cursor];
  }
 
  void popFront() {
    cursor += 1;
  }
 
  @property bool empty() { return cursor >= list.length; }
}



