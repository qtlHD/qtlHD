module qtl.plugins.fourstore.lazycsv;

import std.stdio, std.file;
import std.conv        : to;
import std.c.stdlib    : exit;
import std.array       : split;
import std.string      : chomp, format;

/* Parse a CSV text file into elemidx and rowidx buffers, then use lazy ranges to access the elements, rows and columns  */
struct LazyCsvReader{
  this(string filename, string sep = "\t", string newline = "\n"){           // Index the file, using sep and newlines
    this.filename = filename;
    this.sep      = sep;
    this.newline  = newline;
    filesize      = cast(uint) getSize(filename);
    fp            = new File(filename,"rb");
    if(!exists(filename) || !filename.isFile) throw(new Exception("No such file: '" ~ filename ~ "'"));         // @t Pjotr: Should I: Error Assert Throw ?
    buffer    = new ubyte[](bsize);
    size_t colcnt = 0;                                                       // Current columns counted
    while(fp.rawRead(buffer)){
      foreach(size_t offset, byte b ; buffer){
        if(cast(char)b == sep[0]){
          elemidx ~= bcnt*bsize + (offset+1);                                // Element offsets
          colcnt++;
        }else if(cast(char)b == newline[0]){
          elemidx ~= bcnt*bsize + (offset+1);
          rowidx ~= (elemidx.length-1);                                      // Row start pointers to element offset
          if(colcnt > ncol) ncol = colcnt;                                   // Remember the maximum number of columns
          colcnt = 0;
        }
      }
      bcnt++;
    }
    if(elemidx[($-1)] != filesize) elemidx ~= filesize;                      // Add trailing ends
    if(rowidx[($-1)] != filesize) rowidx ~= (elemidx.length-1);              // Add trailing ends
  }

  // Get an element from the file by element index
  string getElement(size_t l){
    if(l >= elemidx.length) throw(new Exception("No such element: " ~ to!string(l) ~ " (length: " ~ to!string(elemidx.length) ~ ")"));
    fp.seek(elemidx[l]);                                                    // Seek to the location of the element
    ubyte[] buf = new ubyte[](elemidx[l+1]-elemidx[l]);
    fp.rawRead(buf);
    return chomp(cast(string)buf.dup);
  }

  // Get a row by index
  string[] getRow(size_t l){
    if(l >= (rowidx.length-1)) throw(new Exception("No such row: " ~ to!string(l) ~ " (length: " ~ to!string(rowidx.length - 1) ~ ")"));
    size_t start = elemidx[rowidx[l]];
    size_t next = elemidx[rowidx[l+1]];
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
    for(size_t x = 0; x < (rowidx.length-2); x++){                    // 2 dimensions less because we add 0 and  $
      ubyte[] buf;
      if(rowidx[x]+(l+1) < elemidx.length){                           // Next element SHOULD be there
        size_t start    = elemidx[rowidx[x] + l];
        size_t next     = elemidx[rowidx[x] + (l+1)];
        size_t nextrow  = elemidx[rowidx[x+1]];
        if(start < filesize && next-start > 1 && nextrow > start){
          fp.seek(start);
          buf = new ubyte[](next-start-1);
          if(buf) fp.rawRead(buf);
        }
      }
      col[x] = chomp(cast(string)buf.dup);
    }
    return col;
  }

  string toString(){
    return format("'%s': %d Elements, %d Rows used %d Buffers", filename, elemidx.length, rowidx.length, bcnt);
  }

  void close(){ fp.close(); }

  size_t[] elemidx  = [0];   // Public because the ranges need access to the elements length
  size_t[] rowidx   = [0];   // Public because the ranges need access to the row length
  size_t   ncol     =  0;    // Public because the ranges need access to the number of columns
  private:
    string    filename, sep = "\t", newline = "\n";
    size_t    bcnt     = 0, bsize    = 1_048_576;
    size_t    filesize = 0;
    ubyte[]   buffer;
    File*     fp;
}

