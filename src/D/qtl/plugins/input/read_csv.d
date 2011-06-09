/**
 * Plugin for reading simple CSV files
 *
 *   dmd -unittest qtl/plugins/input/readcsv.d ; ./readcsv
 */

module qtl.plugins.input.read_csv;

import qtl.core.primitives;
import std.stdio;

/** 
 * Read a simple CSV file containing marker names, chromosome nrs, position, 
 * phenotype and genotype - such as the listeria.csv file used in R/qtl.
 *
 * The file is parsed once on class instantiation. Elements can be queried.
 */

class ReadSimpleCSV {

  private File f;

  this(in string fn) {
    f = File(fn,"r");
    char[] buf;
    while (f.readln(buf))
      writeln(buf ~ "xx");
    f.close();
  }
}

unittest {
  auto fn = "../../test/data/input/listeria.csv";
  writeln("Reading CSV file" ~ fn);
  auto csv = new ReadSimpleCSV(fn);
}

void main() { }
