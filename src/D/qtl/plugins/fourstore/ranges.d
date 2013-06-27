module qtl.plugins.fourstore.ranges;

import std.stdio, std.file, std.array;
import qtl.plugins.fourstore.lazycsv;

// A lazy forward range of Rows from a CSV file
struct rangeOfRows{
  LazyCsvReader r;
  size_t i = 0;
  string[] fBuffer;

  void popFront(){ i++; }

  @property ref string[] front(){ 
    fBuffer = r.getRow(i);
    return(fBuffer);
  }

  @property bool empty(){ return i==(r.rowidx.length-2); } // We have [0 and filesize]
}
// Helper function to do lazy iteration by row
rangeOfRows byRow(LazyCsvReader r){ return rangeOfRows(r); }

// A lazy forward range of Columns from a CSV file
struct rangeOfColumns{
  LazyCsvReader r;
  string[] fBuffer;
  size_t i = 0;

  void popFront(){ i++; }

  @property ref string[] front(){
    fBuffer = r.getCol(i); 
    return(fBuffer);
  }

  @property bool empty(){ return i==(r.ncol+1); }
}
// Helper function to do lazy iteration by column
rangeOfColumns byColumn(LazyCsvReader r){ return rangeOfColumns(r); }

//  Example should be turned into unit-test
unittest {
  writeln("Unit test " ~ __FILE__, " : plugins.fourstore.ranges");
  string file = "../../test/data/input/hyper.csv";

  LazyCsvReader r = new LazyCsvReader(file, ",");
  writeln(r);  // Print some information

  foreach(col; r.byColumn()){
    writeln("Col: ", col[0..5]);
  }

  foreach(row; r.byRow()){
    writeln("Row: ", row[0..5]);
  }

  r.close();
}

