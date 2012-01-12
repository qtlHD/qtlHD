# qtlHD programming tutorial

This is a gentle (we hope) introduction to data structures and programming
qtlHD.  If you are a C coder, we believe the largest problem has to do with
templates (i.e. generics). The second largest problem is, probably, concepts of
functional programming. As we tend to minimize OOP, we do not think we get into
problems there.

First, we will run through some concepts by reading and writing qtab files. See
also the file
https://github.com/pjotrp/qtlHD/blob/master/doc/input/mapping_inputs.md for a
description of the tabulated qtab file format.

## Reading qtab phenotypes

Start with the unit tests in read_qtab.d (in directory qtl/plugins/qtab).
Currently we read the CSV with read_csv.d. The implementation of read_csv.d
creates a class ReadSimpleCSV wich contains the data as sub-attributes (i.e.
data.phenotypes), which get referenced in that container, a typical OOP
approach.  We should use Tuples instead - which we are doing with more recent
code. The problem is that objects, like ReadSimpleCSV, contain a lot of extra
bagage, which starts to lead its own life. Tuples are clearly data containers,
and nothing else. Tuples are totally appropriate when the goal is to only pass
along data.

So, in write_qtab.d, we read the CSV data in the contained unittest.
Find the function write_phenotype_qtab and you see we pass the data
types, such as individuals and phenotypenames, as explicit input
parameters into the function. That way we can see in one go what the
interface is like. Predictable:

        void write_phenotype_qtab(P)(File f, string descr, in Individuals individuals, in string[] phenotypenames, in P[][] phenotypes)

if we had passed in the ReadSimpleCSV object as an opaque unit, we would not
really know what elements/attributes the function would use, or even modify.

Next, in read_qtab.d, read_qtab reads the record and returns the phenotype qtab
file into a Tuple. Here we access a phenotype matrix

        auto res = read_qtab!(Phenotype!double)(pheno_fn);
        Phenotype!double[][] pheno = res[0];
        // 1st ind, 1st phenotype
        assert(pheno[0][0].value == 118.317);

please check this code. Is it clear what happens? read_qtab uses a
Phenotype!double type. Later we may implement a call with Phenotype!int, for
example. That is the concept of generics/templates. The read_qtab!type function
gets defined, where we set the parameter *type* at compile time.

In the function read_qtab you can see Phenotype!double is used as a
matrix of values (P[][]). The phenotypes get parsed in the data file section

        if (strip(buf) == "# --- Data Phenotype begin") 

where

        auto res = parse_phenotype_qtab(buf);

parses the input string into substrings (returning a Tuple of the
ID, and the values as an array of strings).

The next one is a bit of functional magic

        P[] ps = std.array.array(map!((a) {return set_phenotype!double(a);})(fields));

The function map!(func)(list) takes a list (or array), and applies func to list,
returning a new list. It is similar to the R apply function.

std.array.array() coerces the return type of map back into an array, for further use.

We return the matrix of Phenotypes in a Tuple. That is the container type to use for
a 'bag' of different data structures.

Finally, we test for values with

        // get the first element of the Tuple (there is only one now):
        Phenotype!double[][] pheno = p_res[0];
        // test 1st ind, 1st phenotype
        assert(pheno[0][0].value == 118.317);

The test passes.

From here on you should know how to use the Phenotype!double matrix
for your code.

## Reading qtab symbols

To read the genotype data we need to map the symbols used (A, B, C, H etc) onto
their observed genotypes. After parsing the Listeria CSV write_qtab writes the
symbol table (the example actually uses the prefab F2 observed genotype set,
also defined in genotype.d):

        # --- qtlHD-in-0.1 Symbol Test
        # --- Symbol Genotype begin
        NA - as None
        A AA as 0,0
        B BB as 1,1
        H AB BA as 0,1
        HorB C as 0,1 1,1
        HorA D as 0,0 0,1
        # --- Symbol Genotype end

This needs to be parsed by read_qtab into the in-memory types we use.
There are three basic types (in genotype.d):

1. TrueGenotype contains the two founder alleles.
2. GenotypeCombinator brings multiple names/encodings (i.e. 'A' and 'AA')
   and the true genotypes together (with 'A' is is 0,0). It maintains
   two lists. One for the names/encodings, the other for the matching,
   or observed, true genotypes (see 1.)
3. ObservedGenotypes holds a list of GenotypeCombinators. Which is, in
   effect, the symbol table.

The unittest in read_qtab reads

        unittest {
        The unittest in read_qtab:

        unittest {
          // Symbol and genotype reader
          alias std.path.buildPath buildPath;
          auto dir = to!string(dirName(__FILE__) ~ sep ~ buildPath("..","..","..","..",
        "..","test","data"));
          auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
          writeln("reading ",symbol_fn);
          auto f = File(symbol_fn,"r");
          auto symbols = read_genotype_symbol_qtab(f);
          assert(symbols.decode("A") == symbols.decode("AA"));
          assert(to!string(symbols.decode("A")) == "[(0,0)]");
          assert(to!string(symbols.decode("H")) == "[(0,1)]");
          assert(to!string(symbols.decode("HorA")) == "[(0,0), (0,1)]", to!string(symbo
        ls.decode("HorA")));
        }

which parses the symbol file, and loads above data structures. The three
asserts do a symbol lookup and find the observed genotypes, as stored
in the GenotypeCombinators.

There is no magic here. Parse the genotype table next, but it is really
repeating the above, so every symbol in the table points to its
matching GenotypeCombinator (the list of true genotypes).

It is GenotypeCombinators that we will be using for our algorithms.

## More ...

More to follow...

Copyright (C) 2012 Pjotr Prins <pjotr.prins@thebird.nl> 
