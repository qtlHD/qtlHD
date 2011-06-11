/**
 * Plugin for reading simple CSV files
 *
 *   rm read_csv ; dmd -unittest qtl/plugins/input/read_csv.d qtl/core/primitives.d ; ./read_csv
 */

module qtl.plugins.input.read_csv;

import qtl.core.primitives;
import std.stdio;
import std.regexp;
import std.conv;
import std.string;
import std.path;

/** 
 * Read a simple CSV file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the listeria.csv file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */

class ReadSimpleCSV {

  private File f;
  string phenotypename;
  Marker[] markers;

  this(in string fn) {
    alias std.regexp.split split;

    f = File(fn,"r");
    // read markers
    Marker[] ms;
    foreach (i, mname; split(f.readln(), RegExp("\\s*,\\s*", "i")))
    {
       // writeln(mname);
       Marker m = { name:mname, id:i-1};
       ms ~= m;
    }
    phenotypename = ms[0].name;
    markers = ms[1..$];
    // read chromosome info
    foreach (i, cname; split(f.readln(), RegExp("\\s*,\\s*", "i")))
    {
       auto cname2 = strip(cname);
       // We get a chromosome name, normally a number, or an 'X'.
       // if (cname != "" && cname[0] != 'X' && i > 0) 
         // markers[i-1].chromosome = to!int(cname);
    }
    // read rest
    char[] buf;
    // while (f.readln(buf))
    //  writeln(buf ~ "xx");
    f.close();
  }
}

unittest {
  alias std.path.join join;
  auto fn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","listeria.csv");
  writeln("Unit tests; reading CSV file" ~ fn);
  Marker m2 = { id:2, position:4.8};
  assert(m2.id == 2);
  auto markers = [ m2 ];
  auto data = new ReadSimpleCSV(fn);
  writeln(data.markers.length);
  assert(data.markers.length == 133, to!string(data.markers.length));
  assert(data.phenotypename == "T264");
  writeln(data.markers[0].id);
  assert(data.markers[0].name == "D10M44");
  assert(data.markers[0].id == 0);
  assert(data.markers[1].id == 1);
}

void main() { }
