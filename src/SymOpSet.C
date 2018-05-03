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
// SymOpSet.C:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include "SymOpSet.h"
#include <iomanip>
using namespace std;

SymOpSet::SymOpSet(void) {
  nsym_ = 0;
}

SymOpSet::~SymOpSet() {
}

void SymOpSet::generateOps(char const* name) {

  if (name == "cubic") {

    nsym_ = 48;
    symset_.resize(nsym_);

    symset_[0] = new SymOp(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);       
    symset_[1] = new SymOp(-1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0);     
    symset_[2] = new SymOp(-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,-1.0);     
    symset_[3] = new SymOp(1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,-1.0);     
    symset_[4] = new SymOp(0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0);      
    symset_[5] = new SymOp(0.0,-1.0,0.0,-1.0,0.0,0.0,0.0,0.0,-1.0);    
    symset_[6] = new SymOp(0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0);      
    symset_[7] = new SymOp(0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0);      
    symset_[8] = new SymOp(0.0,0.0,1.0,0.0,-1.0,0.0,1.0,0.0,0.0);      
    symset_[9] = new SymOp(0.0,0.0,-1.0,0.0,-1.0,0.0,-1.0,0.0,0.0);    
    symset_[10] = new SymOp(0.0,0.0,-1.0,0.0,1.0,0.0,1.0,0.0,0.0);      
    symset_[11] = new SymOp(0.0,0.0,1.0,0.0,1.0,0.0,-1.0,0.0,0.0);      
    symset_[12] = new SymOp(-1.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0);      
    symset_[13] = new SymOp(-1.0,0.0,0.0,0.0,0.0,-1.0,0.0,-1.0,0.0);    
    symset_[14] = new SymOp(1.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0);      
    symset_[15] = new SymOp(1.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0);      
    symset_[16] = new SymOp(0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,0.0);       
    symset_[17] = new SymOp(0.0,0.0,-1.0,-1.0,0.0,0.0,0.0,1.0,0.0);     
    symset_[18] = new SymOp(0.0,0.0,-1.0,1.0,0.0,0.0,0.0,-1.0,0.0);     
    symset_[19] = new SymOp(0.0,0.0,1.0,-1.0,0.0,0.0,0.0,-1.0,0.0);     
    symset_[20] = new SymOp(0.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0);       
    symset_[21] = new SymOp(0.0,-1.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0);     
    symset_[22] = new SymOp(0.0,-1.0,0.0,0.0,0.0,1.0,-1.0,0.0,0.0);     
    symset_[23] = new SymOp(0.0,1.0,0.0,0.0,0.0,-1.0,-1.0,0.0,0.0);     
    symset_[24] = new SymOp(-1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,-1.0);    
    symset_[25] = new SymOp(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,-1.0);      
    symset_[26] = new SymOp(1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0);      
    symset_[27] = new SymOp(-1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);      
    symset_[28] = new SymOp(0.0,-1.0,0.0,-1.0,0.0,0.0,0.0,0.0,1.0);     
    symset_[29] = new SymOp(0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0);       
    symset_[30] = new SymOp(0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,-1.0);     
    symset_[31] = new SymOp(0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0);     
    symset_[32] = new SymOp(0.0,0.0,-1.0,0.0,1.0,0.0,-1.0,0.0,0.0);     
    symset_[33] = new SymOp(0.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,0.0);       
    symset_[34] = new SymOp(0.0,0.0,1.0,0.0,-1.0,0.0,-1.0,0.0,0.0);     
    symset_[35] = new SymOp(0.0,0.0,-1.0,0.0,-1.0,0.0,1.0,0.0,0.0);     
    symset_[36] = new SymOp(1.0,0.0,0.0,0.0,0.0,-1.0,0.0,-1.0,0.0);     
    symset_[37] = new SymOp(1.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0);       
    symset_[38] = new SymOp(-1.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0);     
    symset_[39] = new SymOp(-1.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0);     
    symset_[40] = new SymOp(0.0,0.0,-1.0,-1.0,0.0,0.0,0.0,-1.0,0.0);    
    symset_[41] = new SymOp(0.0,0.0,1.0,1.0,0.0,0.0,0.0,-1.0,0.0);      
    symset_[42] = new SymOp(0.0,0.0,1.0,-1.0,0.0,0.0,0.0,1.0,0.0);      
    symset_[43] = new SymOp(0.0,0.0,-1.0,1.0,0.0,0.0,0.0,1.0,0.0);      
    symset_[44] = new SymOp(0.0,-1.0,0.0,0.0,0.0,-1.0,-1.0,0.0,0.0);    
    symset_[45] = new SymOp(0.0,1.0,0.0,0.0,0.0,1.0,-1.0,0.0,0.0);      
    symset_[46] = new SymOp(0.0,1.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0);      
    symset_[47] = new SymOp(0.0,-1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0);
  }
  else {
    cout << "SymOpSet.generateOps, name " << name << " undefined."  << endl;
  }

}

void SymOpSet::convertOpsToXtal(const UnitCell& uc) {
  for (int i=0; i<symset_.size(); i++)
    symset_[i]->convertToXtal(uc);
}

void SymOpSet::printXtal(ostream& os) { 
  for (int i=0; i<symset_.size(); i++)
    symset_[i]->print(os);
}

void SymOpSet::addOp(SymOp* addsym) {
  symset_.push_back(addsym);
  nsym_++;
  assert(nsym_ == symset_.size());
}

void SymOpSet::removeOp(SymOp* delsym) {
  int delindex = -1;
  for (int i=0; i<symset_.size(); i++) {
    if (symset_[i] == delsym)
      delindex = i;
  }
  if (delindex > -1) {
    symset_.erase(symset_.begin()+delindex,symset_.begin()+delindex+1);
    nsym_--;
    assert(nsym_ == symset_.size());
  }
}

void SymOpSet::clear(void) {
  for (int i=0; i<symset_.size(); i++)
    if (!symset_[i] == 0)
      delete symset_[i];
}
