module qtl.plugins.csv.col_range;

import std.stdio, std.file, std.array;
import std.algorithm : find;
import std.string;
import std.traits;
import qtl.plugins.csv.lazy_read_csv;

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

unittest{
  writeln("Unit test " ~ __FILE__, " : plugins.csv.col_range");
  string file = "../../test/data/input/hyper.csv";
  LazyCsvReader r = LazyCsvReader(file, ",");
  writeln(r);  // Print some information
  size_t cnt = 0;
  foreach(col; r.byColumn()){
    if(col.length > 5 && cnt < 5) writeln("Col: ", col[0..5]);
    cnt++;
  }
  r.close();
}
