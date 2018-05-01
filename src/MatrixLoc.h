////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// MatrixLoc.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef MATRIXLOC_H
#define MATRIXLOC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class MatrixLoc : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "matrix_loc"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 3 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> matrix_loc must have two arguments! </ERROR>" << endl;
      return 1;
    }
    
    int mb = atoi(argv[1]);
    int nb = atoi(argv[2]);
    s->wf.set_local_block(mb,nb);
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.cell_mass;
     return st.str();
  }

  MatrixLoc(Sample *sample) : s(sample) {}
};
#endif
