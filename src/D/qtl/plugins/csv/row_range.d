module qtl.plugins.csv.row_range;

import std.stdio, std.file, std.array;
import std.algorithm : find;
import std.string;
import std.traits;
import qtl.plugins.csv.lazy_read_csv;

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

  @property bool empty(){ return i==(r.rowidx.length-1); } // We have [0 and filesize]
}
// Helper function to do lazy iteration by row
rangeOfRows byRow(LazyCsvReader r){ return rangeOfRows(r); }

// Transform a range of string elements into a array of T, header allows to skip certain rows/columns
T[] rangeToTypeArray(T)(string[] elements, size_t[] header = [], string[] missing = ["-","NA"]){
  T[] result;
  foreach(size_t idx, elem; elements){
    if(find(header, idx) == []){                  // Not a header
      if(find(missing, elem) != []){                // Missing
        result ~= T.init;
      }else if(isNumeric!T && isNumeric(elem)){     // Convert to Numeric
        result ~= to!T(elem);
      }else if(isSomeChar!T && !isNumeric(elem)){   // Convert to character
        if(elem.length > 0) result ~= to!T(elem[0]);
      }else if(!isNumeric!T && !isNumeric(elem)){   // Convert to string
        result ~= to!T(elem);
      }else{                                        // Unknown make it missing
        result ~= T.init;
      }
    }
  }
  return result;
}

unittest{
  writeln("Unit test " ~ __FILE__, " : plugins.csv.row_range");
  string file = "../../test/data/input/hyper.csv";

  LazyCsvReader r = LazyCsvReader(file, ",");
  writeln(r);  // Print some information
  foreach(row; r.byRow()){
    if(row.length > 5) writeln("Row: ", row[0..5]);
  }
  r.close();
}

