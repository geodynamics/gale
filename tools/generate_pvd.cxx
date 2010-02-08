/* A simple program to generate field and particle pvd files 

   You will need Boost (www.boost.org) to compile this.  Almost any
   version should do.  On Debian, you compile it with

     g++ generate_pvd.cxx -o generate_pvd -lboost_filesystem

   To generate a static binary, compile it with

     g++ generate_pvd.cxx -o generate_pvd -lboost_system -lboost_filesystem -static

   Usage: generate_pvd NAME TYPE START END STEP

   This will output two files, fields.pvd and particles.pvd.  Put them
   in the same directory as the other fields and particle files.  In
   paraview, open the pvd files.  You will be able use the movie
   controls (play, pause, step forward/backward, jump to
   beginning/end) to examine the time series.

**  Copyright (C) 2009, California Institute of Technology

**  This software is free software; you can redistribute it and/or
**  modify it under the terms of the GNU General Public License as
**  published by the Free Software Foundation; either version 2 of the
**  License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
**  02110-1301 USA
*/


#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <cstdlib>
#include <iostream>

namespace fs=boost::filesystem;
using namespace std;


int main(int argc, char *argv[])
{
  if(argc<5)
    {
      cerr << "Not enough arguments\n"
           << "Usage: generate_pvd [name] [type] [start] [end] [step]\n";
      exit(1);
    }
  string name(argv[1]);
  string suffix(argv[2]);
  if(suffix!="u" && suffix!="s")
    {
      cerr << "Wrong suffix: " << suffix << "\n"
           << "Only 'u' and 's' allowed\n";
      exit(1);
    }

  int start=atoi(argv[3]);
  int end=atoi(argv[4]);
  int step=1;
  if(argc>5)
    step=atoi(argv[5]);
  
  fs::ofstream pvd(name+".pvd");

  pvd << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
      << "  <Collection>\n";
  for(int i=start; i<=end; i+=step)
    {
      pvd << "    <DataSet timestep=\"" << i
          << "\" file=\"" << name << ".";
      pvd.width(5);
      pvd.fill('0');
      pvd << i << ".pvt" << suffix << "\"/>\n";
    }
  pvd << "  </Collection>\n"
      << "</VTKFile>\n";
}
