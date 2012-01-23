# Mapping inputs

Version 0.10 (December 2011)

This document describes the new qtab standard. qtab is a human readable tab
delimited file format for QTL data.

In general, file standards for QTL mapping are a muddle. The basic premises,
however, are pretty straightforward. For qtlHD, we propose a 'final'
tab-delimited format that will also have a binary (and perhaps XML) counterpart
representation. A choice for tab-delimited (or TSV, tab separated value) files
is made to secure human readability (and editability!) of data. 

Since TSV files contain no 'grammar' we describe content and set options
through section header information. In principle we attempt to DRY (do not
repeat yourself) and sparse data (no inclusion of zero values). The idea is to
have a simple, straightforward standard that can be applied to 95% of cases, is
easily readable by humans, and easily parsable by software.

In addition, the new tab delimited format should be versioned, extensible,
handle names with spaces (except tabs), and allow some validation of section
sizes. Section data is either TSV, or JSON. JSON is used for tree like data
structures.  Finally, we support both multi-file and single-file
representations of the same data. The single file representation is simply 
a concatenation of other files.

Note: The qtlHD project is the single authority for the qtab file standard.
qtlHD, itself, acts as a reference implementation.

Each section (or file) is described individually in

* Symbol             - symbols used to describe (observed) genotypes, 
                       markers, founders, individuals and phenotypes
* Founder            - Founder and cross layout
* Genotype           - marker x individual genotypes
* Marker map         - marker (sex) chromosome positions
* Phenotype          - individual x phenotypes (covariates are
                       included as phenotypes)
* Manifest           - optional

Additional files will match one of these formats.  For example, we may
need a "Phenotype covariate" file, containing covariates that describe
the phenotypes (for example, the age for a given body weight
phenotype, or the chromosome and position and sequence of the
transcript for a gene expression phenotype).  This phenotype covariate
file can be in the same format as the phenotype file, but instead of
individual x phenotype, it will be phenotype x "phenotype covariate".

Before describing the formats we describe parser logic, and the general
file/section header layout.

# Parser logic

Any line that starts with a Hash # symbol is considered a comment. Meta-data
for the parser starts with `# --- <Command>` - so hash and three dashes.

A comment may start also later in a line. This is signified by multiple
spaces, a hash and a space. No tabs after an end-of-line comment are allowed. Not all sections
will allow these end-of-line comments.

# Header layout

A header should be easily recognizable, and give a version. Every file 
or file section should start with a commented line (using spaces)

      # --- qtlHD-in-x.x Sectionname Description

where `x.x` is the file format version and `Sectionname` represents the section name. The 
description is a free string. Valid section headers are, for example,

      # --- qtlHD-in-0.9 Genotype Mouse CC JAX
      # --- qtlHD-in-0.9 Phenotype Mouse CC JAX - Our experiment/date

This line is used to recognise and start parsing sections. After the section
header line, a header is recognized by 

      # --- Header begin

and ends with

      # --- Header end

Valid data sections are defined by each format. Typically a data section starts
and ends with

      # --- Data Name begin
      # --- Data Name end

This is an extensible data format. To avoid getting all different formats,
again, try to fit your data model in an existing format before adding a new
one. And note again that qtlHD is the reference implementation.

# The symbol section

The symbol section, or file, substitutes symbols with values. The first line
should be

      # --- qtlHD-in-x.x Symbol Description

A block should start/end with

      # --- Data Type begin/end

A possible founder symbol table

      # --- Symbol Founder begin
      A/J 0
      C57BL/6J 1
      129S1/SvImJ 2
      ...
      # --- Symbol Founder end

When leaving out the numbers, they will be created automatically. So above
is equal to:

      # --- Symbol Founder begin
      A/J
      C57BL/6J
      129S1/SvImJ
      ...
      # --- Symbol Founder end

Likewise, individuals, phenotypes, markers and chromosomes can have symbol to
value substitutions.

Observed genotypes have a somewhat more complex symbol to value substitution.
To translate a genotype symbol to multiple possible (observed) true genotypes,
a possible genotype symbol table:

      # --- Symbol Genotype begin
      NA - as None        
      A as 0,0
      B BB as 1,1
      C CC as 2,2         # aliases
      AB as 1,0           # directional
      BA as 0,1
      AC CA as 0,2 2,0    # not directional
      D as 3              # expands to 3,3
      E as A/J            # Founder symbol expands to 0,0
      ...
      # --- Symbol Genotype end

A symbol is always a string, and can contain spaces (but no tabs). Note that
remarks are also allowed, starting with # and a space.  The `as` keyword splits
aliased symbols and multiple genotypes (`as` is a reserved symbol). 
The genotype numbers refer to the founders. If founders have symbols,
these symbols can be used too. Symbols are simply expanded to their numeric
values.

`None` is a reserved symbol. In data files `NA` is the identifier for a missing
value, or `None`. You can add symbols for `None`. In the above example we allow dash `-` to 
represent a missing value.

Other reserved symbols are `True` and `False` (and, as mentioned
above, `as`).

An example of a symbol reader can be found [here](https://github.com/pjotrp/qtlHD/blob/master/test/data/regression/test_symbol.qtab). Reader and writer are [here](https://github.com/pjotrp/qtlHD/tree/master/src/D/qtl/plugins/qtab). 

# The founder section

A founder section, or file, contains information about the founders useful for
algorithms.  Possibilities are cross description and allele frequencies.

The file, or section, starts with

      # --- qtlHD-in-x.x Founder Description

An example of founder allele frequencies, the probability at each marker
position

        # --- Data Frequencies begin
        # Founder F1 F2 F3 F4
        M1
        A 0.25 0.25 0.25 0.25
        B 0.20 0.30 0.20 0.30
        M2
        A 0.0  0.0  1.0  0.0
        B 0.0  0.0  0.5  0.5
        ...
        # --- Data Frequencies end

*Note: support for the qtab founder section is not yet implemented in qtlHD*

# The genotype section

The genotype section, or file, contains markers (columns) x individuals (rows),
giving the (observed) genotypes. The values can either be numbers, or symbols
(but not both). 

The file, or section, starts with

      # --- qtlHD-in-x.x Genotype Description

Valid genotypes could be (tab delimited):

        # --- Data Observed begin
        #     M1   M2    M3   M4    M5      M6   
        Ind1   A   AorB   1   1,1   0,1   0,1|1,1
        Ind2   H   AorH   NA  1,1   1,1     1,1
        ...
        # --- Data Observed end

where the last marker of Ind1 represents, for example, an AorH. Note the use of
NA - which is the standard name in qtlHD for a missing value. Also note the
treatment of marker names. Marker names are normally listed in the symbol
table, or are deduced from the marker map file. The names in the columns are
normally ignored (as they start with a Hash #). 

The pipe letter `|` in `0,1|1,1` acts as an 'or' combinator for numeric values.
Do not combine symbols in that way. `A|B`, for example, is illegal. For that
case define a new symbol, e.g. `AorB`.

In the header properties can be defined. Current properties are Directional 
(default `True`), which assumes the genotype is directional, e.g. `0,1` differs
from `1,0`.

For example

        # --- Set Genotype begin
        Directional True                  # default
        # --- Set Genotype end

An example of a genotype reader can be found [here](https://github.com/pjotrp/qtlHD/blob/master/test/data/regression/test_genotype.qtab). Reader and writer are [here](https://github.com/pjotrp/qtlHD/tree/master/src/D/qtl/plugins/qtab). 


# The marker map section

The marker map contains a list of markers and their chromosome + locations.
For example,

        # --- Data Location begin
        #   Chr   Pos
        M1   I   100.0
        M2   I   2000.0
        M3   X   100.0
        M4   X   1500.0
        ...
        # --- Data Location end

where, again, marker names and chromosomes can be symbols. Position is in 
cM by default. To override, you can set it to base pairs (bp or Mbp).

        # --- Set Location begin
        Position cM                       # cM (default), bp or Mbp
        # --- Set Location end

*Note: support for the qtab marker map section is not yet implemented in qtlHD*

# The phenotype section

The phenotype section, or file, contains a list of individuals and their phenotypic values.
For example,

        # --- Data Phenotype begin
        #    Sex   P1      P2
        Ind1   X   100.0   45
        Ind2   I   2000.0  46
        ...
        # --- Data Phenotype end

Where the phenotype names are listed in the symbols file. The optional 
name row is ignored by qtlHD as it starts with a Hash symbol.

Phenotypic values can have types. For example, floating point (which is the default
type), integer, binary, discrete, and perhaps ranges. These are defined in the
header, using the phenotype names:

        # --- Type Phenotype begin
        Sex Discrete                        # types are automatically found from data
        P1  Float                           # double precision value with floating point (default)
        P2  Integer
        P3  Ranges 0..1 1..4 4..112 112..$  # valid ranges
        P4  Discrete  T,N,S                 # predefined discrete types
        P5  Integer 1..$                    # only values larger than 1
        P6  Percentage 0..100               # Percentages with range (floating point)
        ...
        # --- Type Phenotype end

Note that types are not symbols. A symbol, defined in the symbol section,
represents a simple value expansion. A type says something about possible
values. When a type is illegal (say a float for an int), or a type falls
outside a predefined set, which can be checked with `P2`, `P3`, `P4`, `P5`, and `P6`, the
software should throw an error.

In addition, every Phenotype section can have a header with batch properties. A
number of properties are standardized, such as Date, Time, Location, Author,
Temperature.  Other fields can be added freely. So

        # --- Set Phenotype begin
        Id batch001
        Date 20111027
        Time 10:11:00
        Location Nema lab
        Author John Smith
        Temperature 25C
        Remark Humid day
        Daylight 8.5 hours
        Week 3
        ...
        # --- Set Phenotype end

Each 'Set' belongs to the 'Data' section. They are tied together in another
table named Property, which mirrors the Data table using the `Id` field:

        # --- Property Phenotypes begin
        #     P1        P2
        Ind1  batch001  batch002
        Ind2  batch002  batch001
        ...
        # --- Data Phenotypes end


This allows phenotype measurements to be split out according to other
parameters. 

An example of a phenotype reader can be found [here](https://github.com/pjotrp/qtlHD/blob/master/test/data/regression/test_phenotype.qtab). Reader and writer are [here](https://github.com/pjotrp/qtlHD/tree/master/src/D/qtl/plugins/qtab). Note that qtab phenotype is only partially implemented in qtlHD.

# Manifest section

The optional Manifest section or file describes the contents and locations
of the other sections and files, including an MD5 checksum. For example:

      # --- qtlHD-in-x.x Manifest Description
      # --- Sections begin
      Symbol symbols.tab 72c66969011ef06999007c4d38ab3c4e
      Genotype data.tab d4ba9bc0e1e72135c61d4e6699bd4197
      Phenotype data.tab d4ba9bc0e1e72135c61d4e6699bd4197
      ...
      # --- Sections end

Note that, in this example, two of the sections are together in one file, and
share the same MD5.

# qtlHD filenaming conventions

The following filename extensions are used by convention

  filename.qtab    - for the textual TSV edition
  filename.qbin    - for the binary edition 
  filename.qbin.gz - for the binary edition, zipped 
  filename.qxml    - for the xml edition

When sections are in separate files, prefix the extension with _sectionname
(lower case). Examples are

  filename_genotype.qtab    - for the textual TSV edition
  filename_phenotype.qtab  
  filename_symbol.qtab  
  filename_founder.qtab  
  filename_location.qtab  
  filename_manifest.qtab  
  filename_genotype.qbin    - for the binary edition 

*Note: support for the qtab manifest section is not yet implemented in qtlHD*

