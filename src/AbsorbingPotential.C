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
// AbsorbingPotential.C
//
// YY: complex absorbing potential (ref. phys. stat. sol. (b) 243, No.5, 1121-1138(2006))
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "AbsorbingPotential.h"

#include "Basis.h"
#include "FourierTransform.h"
#include "D3vector.h"
#include <cassert>
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
AbsorbingPotential::AbsorbingPotential(ChargeDensity& cd):
    cd_(cd), vbasis_(*cd_.vbasis()), vft_(*cd_.vft())
{
   initialize();
}

////////////////////////////////////////////////////////////////////////////////
void AbsorbingPotential::initialize()
{
  R0_ = 7.0;
  deltaR_ = 14.0;
  W0_ = 4.0 / 27.2114;
}

////////////////////////////////////////////////////////////////////////////////
void AbsorbingPotential::update(vector<complex<double>> & vabs)
{
  int np0 = vft_.np0();
  int np1 = vft_.np1();
  int np2 = vft_.np2();
  const complex<double> I(0.0,1.0);
  
  assert(vbasis_.cell().volume() > 0.0);
  UnitCell cell = vbasis_.cell();
  const int np012loc = vft_.np012loc();
  double R0 = R0_;
  double deltaR = deltaR_;
  double W0 = W0_;

  int idx0 = vft_.np0() * vft_.np1() * vft_.np2_first();
  int idxx, i, j, k;
  D3vector r;

  for ( int ir = 0; ir < np012loc; ir++ ) {
    idxx = idx0 + ir;
    k = idxx / ( np0 * np1 );
    idxx = idxx - ( np0 * np1 ) * k;
    j = idxx / np0;
    idxx = idxx - np0 * j;
    i = idxx;
    //r = ( cell.a(0) / np0 * i
    //   + cell.a(1) / np1 * j
    //   + cell.a(2) / np2 * k ) - (cell.a(0) + cell.a(1) + cell.a(2)) / 2.0 ;
    r = ( cell.a(0) / np0 * i
       + cell.a(1) / np1 * j
       + cell.a(2) / np2 * k ) ;
    //cout << r << endl;
    //r_length = length(r);
    cell.fold_in_ws(r);
    if (length(r) < R0 || length(r) > R0+deltaR) {
      vabs[ir] = complex<double>(0.0, 0.0);
    } else 
    {
      // length(r) > R0
      vabs[ir] = -I * W0 * (length(r)-R0)/deltaR;
    }
    
  }
}

