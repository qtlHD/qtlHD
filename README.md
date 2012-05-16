# R/qtlHD

R/qtlHD (or qtlHD) is the next generation implementation of
[R/qtl][rqtl]. qtlHD is aimed at QTL analysis of high-density,
high-throughput data; for example for RNA-seq.

## Features

* qtlHD is written in the high performance [D language](http://dlang.org/). D is a modern 
  strictly typed compiled language with great support for concurrancy. D can be 
  bound against R, Ruby, Python, Perl, etc. D runs on Linux, OS/X and Windows.
* qtlHD currently does scanone
* qtlHD comes with a new flexible data standard for QTL mapping, named qtab. Qtab
  should be able to describe most datasets, see 
  [qtab documentation](https://github.com/kbroman/qtlHD/blob/master/doc/input/qtab.md).

## TODO

Before a full release:

* Create a standalone executable for scanone
* Create an interface for R and [bioruby](https://github.com/pjotrp/bioruby-qtlHD).

## Compiling qtlHD

To compile qtlHD a recent edition of the [D compiler][D] is needed, currently:

    DMD64 D Compiler v2.059

After D installation, try the build and test script (which needs
[rake][rake]):

    ./test_all.sh

## Platforms

qtlHD is developed and tested on Mac OS X, Microsoft Windows and Debian Linux.

## License

qtlHD source code is published under the most liberal open source BSD
License - see the LICENSE file in the [source tree][source].

## Source code

The source code for qtlHD can be found on [github][source].

## Contributors

* Karl W. Broman 
* Pjotr Prins
* Danny Arends

[rqtl]: http://www.rqtl.org/
[D]: http://dlang.org/
[rake]: http://rake.rubyforge.org/
[source]: https://github.com/kbroman/qtlHD
