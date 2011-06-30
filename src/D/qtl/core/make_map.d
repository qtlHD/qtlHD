/**
 * Module for making a new map with inserted inter-marker locations
 */

module qtl.core.make_map;

import std.container;
import std.stdio;
import std.conv;
import std.string;
import std.path;
import std.exception;
import std.algorithm;
alias std.algorithm.find find;

import qtl.core.primitives;
import qtl.core.genotype;

/**
 * The make_fixed_map function creates a new map for a chromosome, and 
 * inserts Pseudo markers at fixed positions.
 *
 * In R/qtl one function existed for fixed_distance_map, 
 * fixed_distance_map_sex and add_marker_if_single.
 *
 * markers:        List of markers (one chromosome)
 * step:           Step size (in cM)
 * off_end:        Max distance (in cM) past the terminal marker
 *
 * (status: under review)
 */

Ms add_stepped_markers_autosome(Ms)(in Ms markers, Position step=1.0, Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);
  // With R/qtl, if step is zero, the purpose was to add a marker if there is
  // only one. Now, use the function add_marker_if_single instead. Adding
  // markers (at off_end) with step=0 is not supported here. That should also
  // be a separate function.
  auto new_markers = markers.sorted();
  if (markers.list.length > 1) {
    auto first_marker = new_markers.list[0];
    auto minpos = first_marker.get_position();
    auto last_marker = new_markers.list[$-1];
    auto maxpos = last_marker.get_position();
    // between minpos and maxpos add markers at fixed positions
    for (auto npos = minpos ; npos < maxpos ; npos += step) {
      // if marker does not exist, add pseudo marker 
      // (note new_markers is sorted)
      auto pm = new PseudoMarker(npos);
      if (countUntil(new_markers.list, pm) == -1)
        new_markers.add(pm);
    }
    // always add one marker beyond each end (no matter
    // the fixed position distance)
    new_markers.add(new PseudoMarker(minpos - off_end));
    new_markers.add(new PseudoMarker(maxpos + off_end));
    new_markers = new_markers.sorted();
    // FIXME remove Pseudo markers too close to other markers
    // (was this in R/qtl?)
  }
  return new_markers;
}

/**
 * NYI FIXME - implementation of sex chromosome will be done later
 */

Ms add_stepped_markers_sex(Ms)(in Ms markers, Position step=1.0, Position off_end=0.0)
{
  /*
  else { # sex-specific map
    if(stepwidth == "variable") {
      if(off.end > 0) {
        tmp <- colnames(map)
        map <- cbind(map[, 1] - off.end, map, map[, ncol(map)] + off.end)
        dimnames(map) <- list(NULL, c("loc000", tmp, "loc999"))
      }
      if(step == 0)
        return(unclass(map))

      ## Determine differences and expansion vector.
      dif <- diff(map[1, ])
      expand <- pmax(1, floor(dif / step))

      ## Create pseudomarker map.
      a <- min(map[1, ]) + cumsum(c(0, rep(dif / expand, expand)))
      b <- min(map[2, ]) + cumsum(c(0, rep(diff(map[2, ]) / expand, expand)))

      namesa <- paste("loc", seq(length(a)), sep = "")
      namesa[cumsum(c(1, expand))] <- dimnames(map)[[2]]
      map <- rbind(a,b)
      dimnames(map) <- list(NULL, namesa)

      return(unclass(map))
    }
    if(stepwidth == "max") {
      if(step==0 && off.end==0) return(unclass(map))
      if(step==0 && off.end>0) {
        if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
          L1 <- L2 <- 1
        }
        else {
          L1 <- diff(range(map[1,]))
          L2 <- diff(range(map[2,]))
        }

        nam <- colnames(map)
        nmap1 <- c(map[1,1]-off.end, map[1,], map[1,ncol(map)]+off.end)
        nmap2 <- c(map[2,1]-off.end*L2/L1, map[2,], map[2,ncol(map)]+off.end*L2/L1)
        map <- rbind(nmap1, nmap2)
        colnames(map) <- c("loc1", nam, "loc2")
        return(unclass(map))
      }

      if(ncol(map)==1) L1 <- L2 <- 1
      else {
        L1 <- diff(range(map[1,]))
        L2 <- diff(range(map[2,]))
      }

      nam <- colnames(map)

      if(off.end > 0) {
        toadd1 <- c(map[1,1] - off.end, map[1,ncol(map)]+off.end)
        toadd2 <- c(map[2,1] + off.end*L2/L1, map[2,ncol(map)]+off.end*L2/L1)

        neword <- order(c(map[1,], toadd1))
        nmap1 <- c(map[1,], toadd1)[neword]
        nmap2 <- c(map[2,], toadd2)[neword]
      }
      else {
        nmap1 <- map[1,]
        nmap2 <- map[2,]
        toadd1 <- toadd2 <- NULL
      }

      d <- diff(nmap1)
      nadd <- ceiling(d/step)-1
      if(sum(nadd) > 0) {
        for(j in 1:(length(nmap1)-1)) {
          if(nadd[j]>0) {
            toadd1 <- c(toadd1, seq(nmap1[j], nmap1[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
            toadd2 <- c(toadd2, seq(nmap2[j], nmap2[j+1], len=nadd[j]+2)[-c(1,nadd[j]+2)])
          }
        }
      }
      newnam <- paste("loc", 1:length(toadd1), sep="")
      
      toadd1 <- sort(toadd1)
      toadd2 <- sort(toadd2)
      neword <- order(c(map[1,], toadd1))
      nmap1 <- c(map[1,], toadd1)[neword]
      nmap2 <- c(map[2,], toadd2)[neword]
      map <- rbind(nmap1, nmap2)
      colnames(map) <- c(nam, newnam)[neword]

      return(unclass(map))
    }

    minloc <- c(min(map[1,]),min(map[2,]))
    map <- unclass(map-minloc)
    markernames <- colnames(map)

    if(step==0 && off.end==0) return(map+minloc)
    else if(step==0 && off.end > 0) {
      map <- map+minloc
      if(ncol(map)==1) { # only one marker; assume equal recomb in sexes
        L1 <- L2 <- 1
      }
      else {
        L1 <- diff(range(map[1,]))
        L2 <- diff(range(map[2,]))
      }

      nam <- colnames(map)
      nmap1 <- c(map[1,1]-off.end, map[1,], map[1,ncol(map)]+off.end)
      nmap2 <- c(map[2,1]-off.end*L2/L1, map[2,], map[2,ncol(map)]+off.end*L2/L1)
      map <- rbind(nmap1, nmap2)
      colnames(map) <- c("loc1", nam, "loc2")
      return(map)
    }
    else if(step>0 && off.end == 0) {

      if(ncol(map)==1) return(map+minloc)

      a <- seq(floor(min(map[1,])),max(map[1,]),
               by = step)
      a <- a[is.na(match(a,map[1,]))]

      if(length(a)==0) return(map+minloc)

      b <- sapply(a,function(x,y,z) {
          ZZ <- min((seq(along=y))[y > x])
          (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1] }, map[1,],map[2,])

      m1 <- c(a,map[1,])
      m2 <- c(b,map[2,])

      names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
      return(rbind(sort(m1),sort(m2))+minloc)
    }
    else {
      a <- seq(floor(min(map[1,])-off.end),ceiling(max(map[1,])+off.end+step),
               by = step)
      a <- a[is.na(match(a,map[1,]))]
      # no more than one point above max(map)+off.end
      z <- (seq(along=a))[a >= max(map[1,])+off.end]
      if(length(z) > 1) a <- a[-z[-1]]

      b <- sapply(a,function(x,y,z,ml) {
        if(x < min(y)) {
          return(min(z) - (min(y)-x)/diff(range(y))*diff(range(z)) - ml)
        }
        else if(x > max(y)) {
          return(max(z) + (x - max(y))/diff(range(y))*diff(range(z)) - ml)
        }
        else {
          ZZ <- min((seq(along=y))[y > x])
          (x-y[ZZ-1])/(y[ZZ]-y[ZZ-1])*(z[ZZ]-z[ZZ-1])+z[ZZ-1]
        }
        }, map[1,],map[2,], minloc[2])
      m1 <- c(a,map[1,])
      m2 <- c(b,map[2,])
      names(m1) <- names(m2) <- c(paste("loc",a,sep=""),markernames)
      return(rbind(sort(m1),sort(m2))+minloc)
    }
  }
}
*/
  return null;
}

/**
 * Add a marker if there is only one. Ms is a list of markers, and 
 * this function works as long as there markers.list and a
 * markers.add function are defined in the Ms type(!)
 */

Ms add_one_if_single_marker(Ms)(in Ms markers, Position step_right=1.0) {
  enforce(step_right>0);
  auto new_markers = new Ms(markers);
  if (markers.list.length == 1) {
    // create a single new marker to the right
    auto marker = new_markers.list[0];
    auto npos = marker.get_position() + step_right;
    auto pm = new PseudoMarker(npos);
    new_markers.add(pm);
  }
  return new_markers;
}

/**
 * Add a range of markers if there is only one in the list
 *
 * off_end:   The map is filled with stepped markers from
 *            marker.position-off_end to marker.position+off_end. 
 */

Ms add_stepped_if_single_marker(Ms)(in Ms markers, Position step=1.0, Position off_end=0.0) {
  enforce(step>0);
  enforce(off_end>=0);
  auto new_markers = new Ms(markers);
  if (markers.list.length == 1) {
    auto marker = new_markers.list[0];
    auto position = marker.get_position();
    for (auto npos = position-off_end ; npos < position+off_end ; npos += step) {
      auto pm = new PseudoMarker(npos);
      new_markers.add(pm);
    }
  }
  return new_markers;
}

unittest {
  writeln("Unit test " ~ __FILE__);
  auto markers = new Markers!(MarkerRef!F2)();
  markers.list ~= new MarkerRef!F2(10.0);
  assert(markers.list.length == 1);
  // call function. Note we don't have to give the type!
  auto new_markers = add_one_if_single_marker(markers,2.0);
  // make sure the original list did not change...
  assert(markers.list.length == 1, "Length is " ~ to!string(markers.list.length));
  // now test the new list
  assert(new_markers.list.length == 2, "Length is " ~ to!string(new_markers.list.length));
  assert(new_markers.list[1].marker.position == 12.0);
  assert(new_markers.list[1].marker.name == "loc12", new_markers.list[1].marker.name);
  // It should work for any list of markers
  auto markers1 = new Markers!(Marker)();
  markers1.list ~= new Marker(10.0);
  assert(markers1.list.length == 1);
  // call function. Note we don't have to give the type!
  auto new_markers1 = add_one_if_single_marker(markers1,2.0);
  // make sure the original list did not change...
  assert(markers1.list.length == 1, "Length is " ~ to!string(markers.list.length));
  // now test the new list
  assert(new_markers1.list.length == 2, "Length is " ~ to!string(new_markers.list.length));
  assert(new_markers1.list[1].position == 12.0);
  assert(new_markers1.list[1].name == "loc12", new_markers1.list[1].name);

  auto new_markers2 = add_stepped_if_single_marker(markers,1.0,5.0);
  assert(new_markers2.list.length == 11, "Length is " ~ to!string(new_markers2.list.length));
  // --- test autosome pseudomarker insertion
  auto markers2 = new Markers!(Marker)();
  // start with three markers
  markers2.list ~= new Marker(10.0);
  markers2.list ~= new Marker(20.0);
  markers2.list ~= new Marker(30.0);
  auto res2 = add_stepped_markers_autosome(markers2,1.0,1.0);
  auto list = res2.list; // new marker list
  // assert there are no duplicates
  assert(list.length == 23, to!string(list.length)); 
  auto uniq_list = uniq!"a.get_position() == b.get_position()"(list);
  auto pos_list = map!"a.get_position()"(uniq_list);
  assert(equal(pos_list,[9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]), to!string(pos_list));
  auto ts = map!"to!double(a)"(pos_list);
  // nudge mar Result into a double[]
  double ds[];
  foreach (u ; pos_list) { ds ~= u; }  // this can be done better, I am sure
  assert(ds[0] == 9);
  assert(ds.length == 23); // tis proof 
}

