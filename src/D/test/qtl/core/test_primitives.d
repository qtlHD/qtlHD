/******************************************************************************
 * Unit tests for primitives 
 */

module test.qtl.core.test_primitives;

import std.stdio;

import qtl.core.primitives;
import qtl.core.phenotype;

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
  Phenotype p1 = { value:-7.809 };
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

