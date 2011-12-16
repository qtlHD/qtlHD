/**
 * Shared library loader
 **/

module qtl.plugins.renv.libload;

//D 2.0 std imports
private import std.loader;
private import std.stdio;
private import std.conv;

/*
 * Gets a function void* from a HXModule and functionname 
 */
protected void* getFunctionThroughVoid(HXModule shared_library, string functionname){
	void* symbol = ExeModule_GetSymbol(shared_library, functionname);
	if (symbol is null) throw new Exception("Failed to load function address " ~ functionname);
	return symbol;
}

/*
 * Loads a single shared library (dll, so, dylib)
 */
protected HXModule load_library(string library_prefix){
  HXModule shared_library = null;
	version (Windows) {
		shared_library = ExeModule_Load(library_prefix ~ ".dll");
	} else version (linux) {
		shared_library = ExeModule_Load("lib" ~ library_prefix ~ ".so");
	} else version (darwin) {
		shared_library = ExeModule_Load(library_prefix ~".dylib");
	}
  if(shared_library is null) throw new Exception("Unable to find shared library: " ~ library_prefix);
  writefln("Loaded shared library: %s",library_prefix);
  return shared_library;
}

/*
 * Adds the operator call to load_function(T)(lib, name)
 */
package struct function_binding(T) {
  bool opCall(HXModule lib, string name) {
    try{
      *fptr = getFunctionThroughVoid(lib, name);
      return true;
    }catch(Exception e){
      writefln("Cannot bind function: %s",name);
      return false;
    }
  }

  private{
    void** fptr;
  }
}

/*
 * Loads a single function (Needs a live reference to the library)
 */
template load_function(T){
  function_binding!(T) load_function(ref T a) {
    function_binding!(T) res;
    res.fptr = cast(void**)&a;
    return res;
  }
}
