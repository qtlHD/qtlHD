# R/qtlHD

R/qtlHD (or qtlHD) is the next generation implementation of
[R/qtl][rqtl]. qtlHD is aimed at QTL analysis of high-density,
high-throughput data; for example for RNA-seq.

## Compiling qtlHD

To compile qtlHD a recent edition of the [D compiler][D] is needed:

    DMD32 D Compiler v2.053

After D installation, try the build and test script (which needs
[rake][rake]):

    ./test_all.sh

## Platforms

qtlHD is developed and tested on OSX, Microsoft Windows and Debian Linux.

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
[D]: http://www.digitalmars.com/d/index.html
[rake]: http://rake.rubyforge.org/
[source]: https://github.com/kbroman/qtlHD
