module qtl.plugins.fourstore.lazycsv;

import std.stdio, std.file;
import std.conv        : to;
import std.c.stdlib    : exit;
import std.array       : split;
import std.string      : chomp, format;

/* Abort with error code, default: -1 */
void error(in string s, int exitcode = -1){
  writeln(s);
  exit(exitcode);
}

class LazyCsvReader{
  this(string fn, string sep = "\t", string newline = "\n"){           // Index the file
    filename = fn;
    filesize  = cast(uint) getSize(filename);
    fp        = new File(filename,"rb");
    if(!exists(filename) || !filename.isFile) error("No such file");
    buffer    = new ubyte[](bsize);
    size_t colcnt = 0;                                           // Current columns counted
    while(fp.rawRead(buffer)){
      foreach(size_t offset, byte b ; buffer){
        if(cast(char)b == sep[0]){
          elemidx ~= bcnt*bsize + (offset+1);                          // Element offsets
          colcnt++;
        }else if(cast(char)b == newline[0]){
          elemidx ~= bcnt*bsize + (offset+1);
          rowidx ~= (elemidx.length-1);                                // Row start pointers to element offset
          if(colcnt > ncol) ncol = colcnt;                             // Remember the maximum number of columns
          colcnt = 0;
        }
      }
      bcnt++;
    }
    if(elemidx[($-1)] != filesize) elemidx ~= filesize;          // Add trailing ends
    if(rowidx[($-1)] != filesize) rowidx ~= (elemidx.length-1);  // Add trailing ends
  }

  ~this(){ close(); }

  // Get an element by index
  string getElement(size_t l){
    fp.seek(elemidx[l]);
    ubyte[] buf = new ubyte[](elemidx[l+1]-elemidx[l]);
    fp.rawRead(buf);
    return chomp(cast(string)buf.dup);
  }

  // Get a row by index
  string[] getRow(size_t l, string sep = "\t"){
    size_t start = elemidx[rowidx[l]];
    uint next = elemidx[rowidx[l+1]];
    fp.seek(start);
    ubyte[] buf;
    if(start != filesize && next-start > 0){
      buf = new ubyte[next-start];
      fp.rawRead(buf);
    }
    string[] row = split(chomp(cast(string)buf.dup), sep);
    return row;
  }

  // Get a column by index
  string[] getCol(size_t l){
    string[] col = new string[](rowidx.length-2);
    for(size_t x = 0; x < (rowidx.length-2); x++){
      size_t start  = elemidx[rowidx[x]+l];
      uint next     = elemidx[rowidx[x] + (l+1)];
      uint nextrow  = elemidx[rowidx[x+1]];
      ubyte[] buf;
      if(start != filesize && next-start > 1 && nextrow > start){
        fp.seek(start);
        buf = new ubyte[](next-start-1);
        if(buf) fp.rawRead(buf);
      }
      col[x] = chomp(cast(string)buf.dup);
    }
    return col;
  }

  override string toString(){
    return format("'%s': %d Elements, %d Rows used %d Buffers", filename, elemidx.length, rowidx.length, bcnt);
  }

  void close(){ fp.close(); }

  size_t[] elemidx  = [0];
  size_t[] rowidx   = [0];
  size_t   ncol     = 0;
  private:
    string    filename;
    size_t    bcnt     = 0, bsize    = 1_048_576;
    size_t    filesize = 0;
    ubyte[]   buffer;
    File*     fp;
}

/*  Example should be turned into unit-test
void main(string[] args){
  string file = "test.file";
  if(args.length > 1) file = args[1];
  LazyCsvReader r = new LazyCsvReader(file);

  foreach(col; r.byColumn()){
    writeln("Col", col);
  }

  foreach(row; r.byRow()){
    writeln("Row", row);
  }

  writeln(r);  // Print some information
  r.close();
}
*/
