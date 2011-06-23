/**
 * \file read_binary.d - Plugin for writing XBIN files
 *
 * Copyright (c) 2011 Danny Arends
 * Part of the qtlHD package
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 *
 * - Mac / Linux Unittest:
 * dmd -unittest qtl/plugins/input/read_binary.d qtl/core/*.d  
 *
 * - Win32 (No circular dependancies, or multiple definitions of main allowed):
 *
 * dmd -run cdc.d -lib qtl/core/ -ofCore.lib
 * dmd -run cdc.d -unittest qtl/plugins/input/read_binary.d Core.lib
 *
 *
 **/
 
 module qtl.plugins.input.read_binary;

import qtl.core.primitives;
import qtl.core.chromosome;
import qtl.core.phenotype;
import qtl.core.genotype;

import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.file;

import qtl.plugins.input.binary_types;

class XbinReader(XType){
  private File f;
  int nindividuals = 0;
  int nphenotypes = 0;
  int nmarkers = 0;
  
  bool checkBuffer(ubyte[] buffer){
    ubyte[] startprint = buffer[0..2];
    ubyte[] endprint = buffer[buffer.length-2..$];
    if(startprint == b_footprint && endprint == b_footprint){
      return true;
    }
    return false;
  }
  
  int byteToInt(ubyte[] bits, bool little_endian = true ){
    int result = 0;
    if(little_endian){
      for(int n = bits.length-1; n >= 0; n--)
        result = (result << 8) +bits[ n ];
    }else{
      for(int n = 0; n < bits.length; n++)
        result = (result << 8) +bits[ n ];
    }
    return result;
  }
  
  T[] toType(T)(ubyte[] buffer){
    T[] returnbuffer;
    foreach(int i, byte b ; buffer){
      returnbuffer ~= to!T(b);
    }
    return returnbuffer;
  }
  
  int[] getVersion(ubyte[] buffer){
    int[] fileversion = toType!int(buffer[2..5]);
    return fileversion;
  }
  
  int getNumberOfMatrices(ubyte[] buffer){
    int nmatrix = byteToInt(buffer[5..9]);
    return nmatrix;
  }
  
  int getMatrix(ubyte[] buffer, int matrix, int skip){
    int type = byteToInt(buffer[skip..skip+4]);
    int nrow = byteToInt(buffer[(skip+4)..(skip+8)]);
    int ncol = byteToInt(buffer[(skip+8)..(skip+12)]);
    int elem = byteToInt(buffer[(skip+12)..(skip+16)]);
    writefln("matrix %d %d %d %d", matrix, type, nrow, ncol);
    return 100;
  }
  
  this(in string filename){
    assert(getSize(filename) < uint.max);
    ubyte[] inputbuffer = new ubyte[cast(uint)getSize(filename)];
    auto f = new File(filename,"rb");

    f.rawRead(inputbuffer);
    //Header
    writeln("Read: " ~ to!string(getSize(filename)) ~ " bytes");
    writeln("File OK? " ~ to!string(checkBuffer(inputbuffer)));
    writeln("Version: " ~ to!string(getVersion(inputbuffer)));
    writeln("Matrices: " ~ to!string(getNumberOfMatrices(inputbuffer)));
    //Loop through the matrices
    int skip = 9;
    for(int m=0; m<getNumberOfMatrices(inputbuffer); m++){
      skip += getMatrix(inputbuffer, m, skip);
    }
  }
}

unittest{
  writeln("Unit test " ~ __FILE__);
  alias std.path.join join;
  auto infn = dirname(__FILE__) ~ sep ~ join("..","..","..","..","..","test","data","input","multitrait.xbin");
  writeln("  - reading XBIN " ~ infn);
  auto data = new XbinReader!RIL(infn);
}

void main() { }
