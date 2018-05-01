////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
#include <config.h>
#include <vector>
#include "jacobi_eigenvalue.h"
int jade_complex(int maxsweep, double tol, std::vector<ComplexMatrix*> a,
         ComplexMatrix& u, std::vector<std::vector<std::complex<double> > > &adiag);
void jacobi_eigenvalue( int n, double a[], int it_max, double v[],
  double d[], int &it_num, int &rot_num );
