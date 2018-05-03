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
// TDMLWFTransform.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef TDMLWFTRANSFORM_H
#define TDMLWFTRANSFORM_H

#include <vector>
#include <complex>
class SlaterDet;
class UnitCell;
class DoubleMatrix;
class ComplexMatrix; //DCY
#include "D3vector.h"
#include "BasisMapping.h"
#include "UserInterface.h" //DCY
#include <stdlib.h>
#include <string>
#include "Sample.h"
#include "Matrix.h"

class TDMLWFTransform
{
  private:

  const SlaterDet& sd_;
  const UnitCell& cell_;
  const Context& ctxt_;

  BasisMapping bm_;
  //DoubleMatrix* aprox; //proxy double matrix for a's
  std::vector<ComplexMatrix*> a_;  //DCY cosine and sine matrices
  ComplexMatrix* u_;               //DCY orthogonal transformation
  std::vector<std::vector<std::complex<double> > > adiag_; //DCY, diagonal elements
  std::vector<std::vector<double> > adiag_real_;
  std::vector<std::vector<double> > adiag_imag_; 

  SlaterDet *sdcosx_, *sdsinx_,
            *sdcosy_, *sdsiny_,
            *sdcosz_, *sdsinz_;

  public:
  ComplexMatrix* a(int k) { return a_[k]; };

  SlaterDet* sdcosx(void) { return sdcosx_; };
  SlaterDet* sdcosy(void) { return sdcosy_; };
  SlaterDet* sdcosz(void) { return sdcosz_; };
  SlaterDet* sdsinx(void) { return sdsinx_; };
  SlaterDet* sdsiny(void) { return sdsiny_; };
  SlaterDet* sdsinz(void) { return sdsinz_; };

  std::complex<double> adiag(int k, int i) { return adiag_[k][i]; }

  void update(void); //compute matrices for Berry phase and MLWF
  void compute_transform(void);
  void compute_sincos(const int n, const std::complex<double>* f,
    std::complex<double>* fc, std::complex<double>* fs);
  void apply_transform(SlaterDet& sd);

  double spread2(int i, int j);
  double spread2(int i);
  double spread2(void);
  double spread(int i);
  double spread(void);
  D3vector center(int i);
  D3vector dipole(void);

  TDMLWFTransform(const SlaterDet& sd);
  ~TDMLWFTransform(void);
};
#endif
