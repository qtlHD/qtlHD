# qtlHD programming tutorial

This is a somewhat gentle (we hope) introduction to data structures and
programming qtlHD.  If you are a C coder, we believe the largest problem is
templates (aka generics). The second largest problem is, probably, concepts of
functional programming. As we tend to minimize OOP, we do not think we get into
problems there.

First, we will run through some concepts by reading and writing qtab files. See
also the 
[qtab documentation](https://github.com/pjotrp/qtlHD/blob/master/doc/input/qtab.md) for a
description of the tabulated qtab file format, and the implemented qtab [reader and writer](https://github.com/pjotrp/qtlHD/tree/master/src/D/qtl/plugins/qtab).

## Reading qtab phenotypes

Start with the unit tests in [read_qtab.d](https://github.com/pjotrp/qtlHD/blob/master/src/D/qtl/plugins/qtab/read_qtab.d) (in source code directory qtl/plugins/qtab).
Currently we read the CSV with *read_csv.d*. The implementation of read_csv.d
creates a class ReadSimpleCSV which contains the data as sub-attributes (e.g.,
data.phenotypes), which get referenced in that container, a typical OOP
approach.  We should use Tuples instead - which we are doing with more recent
code. The reason to use Tuples is that objects, such as ReadSimpleCSV, contain a lot of extra
bagage, which starts to lead its own life. Tuples are clearly and solely data containers.
Nothing else. Tuples are totally appropriate when the goal is to only pass
data around.

So, for write_qtab.d, we read the CSV data in the contained unittest.
Find the function write_phenotype_qtab and you see we pass the data
types, such as individuals and phenotypenames, as explicit input
parameters into the function. That way we can see in one go what the
interface is like. This interface is predictable:

```D
void write_phenotype_qtab(P)(File f, string descr, in Individuals individuals, in string[] phenotypenames, in P[][] phenotypes)
```

i.e. what you see is what you get.
If we had passed in the ReadSimpleCSV object as an opaque unit, in typical OOP fashion, we would not
really know what elements/attributes the function would use, or even modify.

Next, in read_qtab.d, read_qtab reads the record and returns the phenotype qtab
file into a Tuple. Here we access a phenotype matrix

```D
auto res = read_qtab!(Phenotype!double)(pheno_fn);
Phenotype!double[][] pheno = res[0];
// 1st ind, 1st phenotype
assert(pheno[0][0].value == 118.317);
```

please check this code. Is it clear what happens? read_qtab uses a
Phenotype!double type (I grant the exclamation mark is a bit strange syntax).
Later we may implement a call with Phenotype!int, for example. That is the
concept of generics/templates. The read_qtab!type function gets defined, where
we set the parameter *type* at compile time. The exclamation type combination
makes the implementation generic for a type. This is called 'generics', where
Phenotype is a template, and double is the type. JAVA, Scala, C++ all have
these types of templates. The only difference is that D uses an exclamation
mark to pass the template parameter(s).

In the function read_qtab you can see Phenotype!double is used as a
matrix of values (`P[][]`). The phenotypes get parsed in the data file section

```D
if (strip(buf) == "# --- Data Phenotype begin") 
```

where

```D
auto res = parse_phenotype_qtab(buf);
```

parses the input string into substrings (returning a Tuple of the
ID, and the values as an array of strings).

The next one is a bit of functional magic

```D
P[] ps = std.array.array(map!((a) {return set_phenotype!double(a);})(fields));
```

The function map!(func)(list) takes a list (or array), and applies func to list,
returning a new list. It is similar to the R apply function.

std.array.array() coerces the return type of map back into an array, for further use.

We return the matrix of Phenotypes in a Tuple. That is the container type to use for
a 'bag' of different data structures.

Finally, we test for values with

```D
// get the first element of the Tuple (there is only one now):
Phenotype!double[][] pheno = p_res[0];
// test 1st ind, 1st phenotype
assert(pheno[0][0].value == 118.317);
```

The test passes.

From here on you should know how to use the Phenotype!double matrix
for your code.

## Reading qtab symbols

To read the genotype data we need to map the symbols used (A, B, C, H etc) onto
their observed genotypes. After parsing the Listeria CSV write_qtab writes the
symbol table (the example actually uses the prefab F2 observed genotype set,
also defined in genotype.d):

    # --- qtlHD-in-0.1 Symbol Test
    # --- Genotype Symbol begin
    NA - as None
    A AA as 0,0
    B BB as 1,1
    H AB BA as 0,1
    HorB C as 0,1 1,0 1,1
    HorA D as 0,0 0,1 1,0
    # --- Genotype Symbol end

This needs to be parsed by read_qtab into the in-memory types we use.
There are three basic types (in genotype.d):

1. TrueGenotype contains the two founder alleles.
2. GenotypeCombinator brings multiple names/encodings (e.g., 'A' and 'AA')
   and the true genotypes together (with 'A' is is 0,0). It maintains
   two lists. One for the names/encodings, the other for the matching,
   or observed, true genotypes (see 1.)
3. ObservedGenotypes holds a list of GenotypeCombinators. Which is, in
   effect, the symbol table.

The unittest in read_qtab reads

```D
unittest {
The unittest in read_qtab:

unittest {
  // Symbol and genotype reader
  alias std.path.buildPath buildPath;
  auto dir = to!string(dirName(__FILE__) ~ dirSeparator ~ buildPath("..","..","..","..",
"..","test","data"));
  auto symbol_fn = to!string(buildPath(dir,"regression","test_symbol.qtab"));
  writeln("reading ",symbol_fn);
  auto f = File(symbol_fn,"r");
  auto symbols = read_genotype_symbol_qtab(f);
  assert(symbols.decode("A") == symbols.decode("AA"));
  assert(to!string(symbols.decode("A")) == "[(0,0)]");
  assert(to!string(symbols.decode("H")) == "[(0,1)]");
  assert(to!string(symbols.decode("HorA")) == "[(0,0), (0,1), (1,0)]");
  // It is also possible to encode directly with
  assert(to!string(symbols.decode("0,0")) == "[(0,0)]");
  assert(to!string(symbols.decode("0,1")) == "[(0,1)]");
  assert(to!string(symbols.decode("0,0|0,1")) == "[(0,0), (0,1)]");
}
```

which parses the symbol file, and loads above data structures. The three
asserts do a symbol lookup and find the observed genotypes, as stored
in the GenotypeCombinators.

There is no magic here. Parse the genotype table next, but it is really
repeating the above, so every symbol in the table points to its
matching GenotypeCombinator (the list of true genotypes).

It is GenotypeCombinators that we will be using for our algorithms.

## Reading qtab genotypes

Now we have the genotype symbols we can parse the qtab genotype table. As
described in the qtab format (see link in introduction) this is a matrix of
marker by individual, containing the observed genotypes, e.g. A, B, H, HorB,
etc. The symbols that can be used are defined in the symbol table, which
relates them to the actual genotype. So A would be 0,0 (the genotype of father
and mother is 0), B is 1,1, H is 1,0 (or 0,1), HorA is an observed genotype,
which can be either H or B at the marker location. Obviously you can define
your own symbols. Often C is chosen for HorB, and D for HorA. Any symbol (a
string of text) can be mapped against a list of allowed genotypes. In the
future we may even support polyploidity.

In the previous section we read the symbol table into a ObservedGenotypes
container (defined in genotype.d) named symbols using

```D
auto f = File(symbol_fn,"r");
auto symbols = read_genotype_symbol_qtab(f);
```

Now we read the genotype matrix with

```D
auto f = File(genotype_fn,"r");
auto ret = read_genotype_qtab(f1, symbols);
auto individuals = ret[0];
auto genotype_matrix = ret[1];
// Show the first individual and genotypes
writeln(individuals.list[0].name,genotype_matrix[0]);
```

which returns a Tuple with individuals and accompanying genotypes, in a matrix of references to symbols, i.e. GenotypeCombinators(!). For convenience we named these
Gref. So a matrix of references to (observed) genotype combinators is defined as

```D
Gref genotype_matrix[][]
```

This implies that the matrix contains only references (pointers) to a limited
number of symbols. Which is efficient. In the case of F2 it would be A, B, H,
and perhaps a few combinations theref, such as HorA and HorB. The genotype
matrix simply references these symbols. This allows us to get at the true
genotypes by querying the relevant symbol or ObservedGenotypes container. To 
use this genotype, we can match the decoder:

```D
assert(genotype_matrix[0][0] == symbols.decode("B"));
assert(genotype_matrix[0][3] == symbols.decode("H"));
```

or we can ask for the true genotypes using a symbol

```D
assert(genotype_matrix[0][0] == symbols.decode("B"));
assert(genotype_matrix[0][3] == symbols.decode("H"));
```

or we can query the founder parents directly

```D
assert(genotype_matrix[0][0].list[0].homozygous == true);
assert(genotype_matrix[0][0].list[0].founders[0] == 1);
assert(genotype_matrix[0][0].list[0].founders[1] == 1);
assert(genotype_matrix[0][3].list[0].heterozygous == true);
assert(genotype_matrix[0][3].list[0].founders[0] == 0);
assert(genotype_matrix[0][3].list[0].founders[1] == 1);
```

(we may add some syntactic sugar later). Note both B and H have
one set of parents. If there are more, as in the case of HorB for marker
M you get

```D
assert(genotype_matrix[0][M].list[0].heterozygous == true);
assert(genotype_matrix[0][M].list[0].founders[0] == 0);
assert(genotype_matrix[0][M].list[0].founders[1] == 1);
assert(genotype_matrix[0][M].list[1].homozygous == true);
assert(genotype_matrix[0][M].list[1].founders[0] == 1);
assert(genotype_matrix[0][M].list[1].founders[1] == 1);
```

in other words, we can digest the founder combination for every symbol, which
is the information we gave qtlHD with the symbol table:

    # --- qtlHD-in-0.1 Symbol Test
    # --- Genotype Symbol begin
    NA - as None
    A AA as 0,0
    B BB as 1,1
    H AB BA as 0,1
    HorB C as 0,1 1,0 1,1
    HorA D as 0,0 0,1 1,0
    # --- Genotype Symbol end

At this stage you should know how to use the phenotype and genotype matrices.

## Reading the marker map

The final piece of the QTL mapping puzzle is the marker map. An example can be found
[here](https://github.com/pjotrp/qtlHD/blob/master/test/data/regression/test_marker_map.qtab). This
a marker consists of a name, a chromosome name, and a position on that chromosome. A 
marker info is defined in primitives.d (unsurprisingly) as 

```D
mixin template MarkerInfo() {
  mixin Identity;
  mixin Attributes;
  Chromosome chromosome;      /// Reference to Chromosome
  Position position;          /// Marker position - content depends on map
}
```

where identity gives it a name. Attributes is perhaps the strange one - we use
attributes here as a container for non-generic properties. You can imagine a
marker to be fictional, or some other feature - this we handle through its
attributes. Attributes contain information that does not apply to all
algorithms, but needs to be tracked in some algorithms. Anyway, we can ignore
it for now.

First we parse the qtab marker map in *read_qtab.d*, much the same as to what we did
earlier. Pass in the qtab file (or section) and return a list of markers.

```D
auto markers = read_marker_map_qtab!(Marker)(marker_map_fn);
assert(markers[0].name == "D10M44");
assert(markers[0].chromosome.name == "1");
assert(markers[3].position == 40.4136);
```

The file parser looks like this

```D
auto read_marker_map_qtab(M)(string fn) {
  M ret_ms[];
  auto f = File(fn,"r");
  scope(exit) f.close(); // always close the file on function exit
  string buf;
  while (f.readln(buf)) {
    if (strip(buf) == "# --- Data Location begin") {
      while (f.readln(buf)) { 
        if (strip(buf) == "# --- Data Location end")
           break;
        if (buf[0] == '#') continue;
        auto res = parse_marker_qtab(buf);
        auto name = res[0];
        auto cname = res[1];
        auto c = new Chromosome(cname);
        auto pos = res[2];
        auto marker = new Marker(c,pos,name);
        ret_ms ~= marker;
      }
    }
  }
  return ret_ms;
}
```

The line parser reads 

```D
Tuple!(string, string, double) parse_marker_qtab(string line) {
  auto fields1 = split(line,"\t");
  auto fields = std.array.array(map!"strip(a)"(fields1));  // <- note conversion to array
  auto name = (fields.length > 0 ? strip(fields[0]) : null);
  auto chromosome = (fields.length > 1 ? fields[1] : null);
  auto position = (fields.length > 2 ? to!double(strip(fields[2])) : MARKER_POSITION_UNKNOWN);
  return tuple(name,chromosome,position);
}
```

Pretty straightforward again. Note, however, that the line parser does very little 
syntax checking. Also a new chromosome object gets created with evey marker object - 
maybe not what we want. Also, the allocation of 

```D
M ret_ms[];
```

is dynamically allocated, and may slow things down with appending markers to
really large sets. Anyway, no worries for now.

Finally, for the purists, in the current implementation we see some repetition in
the parsers. It maybe we reduce the code repetition later (DRY).

More to follow...

Copyright (C) 2012 Pjotr Prins <pjotr.prins@thebird.nl> 
