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
  
  bool checkBuffer(byte[] buffer){
    byte[] startprint = buffer[0..2];
    byte[] endprint = buffer[buffer.length-2..$];
    if(startprint == b_footprint && endprint == b_footprint){
      return true;
    }
    return false;
  }
  
  T[] toType(T)(in byte[] buffer){
    T[] returnbuffer;
    foreach(int i, byte b ; buffer){
      returnbuffer ~= to!T(b);
    }
    return returnbuffer;
  }
  
  int[] getVersion(byte[] buffer){
    int[] fileversion = toType!int(buffer[2..5]);
    return fileversion;
  }
  
  this(in string filename){
    assert(getSize(filename) < uint.max);
    byte[] inputbuffer = new byte[cast(uint)getSize(filename)];
    auto f = new File(filename,"rb");
    int byte_cnt=0;
    int buffer_cnt=0;
    f.rawRead(inputbuffer);
    writeln("Read: " ~ to!string(getSize(filename)) ~ " bytes");
    writeln("File OK? " ~ to!string(checkBuffer(inputbuffer)));
    writeln("Version: " ~ to!string(getVersion(inputbuffer)));
    
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
