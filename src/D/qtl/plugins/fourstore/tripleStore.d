module qtl.plugins.fourstore.TripleStore;

import std.stdio, std.math, std.conv, std.socket, std.string;

bool has(T,U)(T[U] array, U key){
  foreach(U k; array.byKey()){
    if(k==key){ return(true); }
  }
  return false;
}

/** Triple store class object */
class TripleStore{
  public:
    /** Constructor for the TripleStore object */
    this(string h = "localhost", ushort p = 9000){ 
      host = h; port = p;
    }

    /** Query the triple store with query string query, supported formats: text, json, sparql */
    string query(string query, string format="text", uint buffersize = 1024){
      query = addPrefixes(query);
      write("GET "~url~"?query="~query~"&output="~format~" "~protocol~"\r\nHost: "~host~"\r\n\r\n");
      return readOutput(buffersize);
    }

    /** Connect to a triple store */
    bool connect(){
      try{
        handle = new TcpSocket();
        handle.connect(new InternetAddress(host, cast(int)port));
        handle.blocking(true);
      }catch(SocketException ex){
        writefln("Failed to connect to server (%s:%d)", host, port);
        handle = null;
        return false;
      }
      return true;
    }

    /** Disconnect from a triple store */
    void disconnect(){
      if(isAlive()){ handle.close; delete handle; }
    }

    void addPrefix(string prefix, string uri){
      prefixes[prefix] = uri;
    }

    bool isPrefix(string prefix){
      if(prefix[($-1)] == ':'){
        return(has(prefixes, prefix));
      }
      return false;
    }

  private:

    string addPrefixes(string query){
      string prequery;
      foreach(string p; prefixes.byKey()){
        prequery ~= "PREFIX " ~ p ~ " " ~ prefixes[p] ~ "\n";
      }
      return (prequery ~ query);
    }

    bool write(string msg){ // Write a message to server
      if(!isAlive()) return false;
      auto ret = handle.send(msg);
      return (ret < 0);
    }

    bool isAlive(){ // Is our handle still alive
      if(handle !is null) return handle.isAlive;
      return false;
    }

    string readOutput(uint buffersize = 1024){ // Read the output from the server  
      if(!isAlive()) return null;
      string response;
      char[] buf = new char[buffersize];
      int ret;
      while((ret = handle.receive(buf)) > 0){
        response ~= to!string(buf[0 .. ret].dup);
      }
      auto header = response.indexOf("\r\n\r\n"); // Remove the html header
      return response[(header+4) .. ($-1)];       // And the trailing newline
    }

    string[string] prefixes;
    string host     = "localhost"; // Host to connect to
    string url      = "/sparql/";  // Url prefix on host for sparql queries
    string protocol = "HTTP/1.0";  // HTTP version we want to use
    ushort port     = 9000;        // Port to connect on
    Socket handle   = null;        // Socket to server (called handle)
}

void main(string[] args){
  TripleStore store = new TripleStore();
  store.addPrefix(":","<http://www.rqtl.org/ns/#>");
  store.addPrefix("individual:","<http://www.rqtl.org/ns/individual#>");
  store.addPrefix("marker:","<http://www.rqtl.org/ns/marker#>");
  store.addPrefix("location:","<http://www.rqtl.org/ns/location#>");
  store.addPrefix("phenotype:","<http://www.rqtl.org/ns/phenotype#>");

  store.connect();
  writefln("--------------------------------------");
  writefln("%s", store.query("SELECT * WHERE {?s location:chr ?o} LIMIT 50", "text"));
  writefln("--------------------------------------");
  store.disconnect();
}

