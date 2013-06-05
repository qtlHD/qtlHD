# qtlHD data structures

In this document we discuss the data structures internal to qtlHD. I.e.
these structures are foremost available to the D modules. These data structures
mostly mirror the input data formats discussed in doc/input/qtab.md.

## Introduction

The basic design principle is KISS. Ideally. Unfortunately, abstractions and a
strong type system quickly does away with simplicity.

Still we want the algorithms to be as clean as possible. This means we can use
the same algorith with a simple type (say, represent an individual just as a
string containing a name), or a more elaborate type, i.e. an object containing
more data. As long as to!string is supported by both, most use cases can be
covered in a simple fashion.

