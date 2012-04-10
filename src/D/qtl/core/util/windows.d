/**
 * Windows related functions
 * Many thanks to: https://github.com/kyllingstad/ltk/
 */

module qtl.core.util.windows;

version(Windows){
  import core.sys.windows.windows;
  import std.utf: toUTF16z;
  import std.windows.syserror;
  import std.stdio;
  import std.conv;

  extern(Windows){
    LPTCH GetEnvironmentStrings();
    DWORD GetEnvironmentVariableW(LPCWSTR lpName, LPWSTR lpBuffer, DWORD nSize);
    BOOL  SetEnvironmentVariableW(LPCWSTR lpName, LPCWSTR lpValue);
  }
  
 /**
  * Check if the variable 'name' is defined in the environment
  */
  private bool envExists(LPCWSTR name){
    return GetEnvironmentVariableW(name, null, 0) != 0;
  }
  
 /**
  * Set a new value for variable 'name' in the environment
  */
  void setEnv(string name, string value, bool overwrite){
    auto name16z = toUTF16z(name);
    if (!overwrite && envExists(name16z)) return;
    SetEnvironmentVariableW(name16z, toUTF16z(value));
    sysErrorString(GetLastError());
  }
 
 /**
  * Get the string value for variable 'name' in the environment
  */
  string getEnv(string name){
    wchar buffer[1024];
    const value16z = GetEnvironmentVariableW(toUTF16z(name), buffer.ptr, 1024);
    return to!string(buffer[0 .. value16z].idup);
  }
  
 /**
  * Delete the variable 'name' in the environment
  */  
  void unsetEnv(string name){
    auto name16z = toUTF16z(name);
    if (envExists(name16z)){
      SetEnvironmentVariableW(name16z, null);
      sysErrorString(GetLastError());
    }
  }
  
  unittest {
    import std.stdio;
  
    setEnv("test", "aaa", true);
    assert (getEnv("test") == "aaa");

    setEnv("test", "bbb", true);
    assert (getEnv("test") == "bbb");

    setEnv("test", "ccc", false);
    assert (getEnv("test") == "bbb");

    unsetEnv("test");
    assert (getEnv("test") == null);
  }
  
}

