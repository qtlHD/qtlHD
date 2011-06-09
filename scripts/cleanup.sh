#! /bin/sh

if [ ! -d scripts ]; then
  echo "Run from root of repo"
fi
find -name *.o -exec rm -v \{\} \;
find -name *.deps -exec rm -v \{\} \;
