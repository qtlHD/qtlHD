/**
 * Program for converting CSV/CSVr file to XGAP and XGAP binary (xbin)
 */

module qtl.util.csv2xgap;

void main(string[] args){
  if(args.length != 3){
    writeln("Usage: convert in.csvr out.xbin");
  }else{
    writeln("reading CSVR (" ~ args[1] ~ ") to XBIN (" ~ args[2] ~ ")");
    auto data = new CsvrReader!RIL(args[1]);
    auto result = new BinaryWriter!(CsvrReader!RIL,RIL)(data,args[2]);
    writefln("Reduced from: %.2f Kb to %.2f Kb", toKb(args[1]), toKb(args[2]));
  }
}

