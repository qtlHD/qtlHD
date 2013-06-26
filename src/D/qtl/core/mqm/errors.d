module qtl.core.mqm.errors;

import std.c.stdlib : exit;
import std.stdio : stderr, writeln, writefln;

/* Write an warning string to stdout */
void warning(in string s){
  writeln();
  writefln("-Warning: %s\n", s);
}

/* Write an error string to stderr */
void error(in string s){
  stderr.writeln();
  stderr.writefln("-Error: %s\n", s);
}

/* Abort with error code, default: -1 */
void abort(in string s, int exitcode = -1){
  error(s);
  exit(exitcode);
}
