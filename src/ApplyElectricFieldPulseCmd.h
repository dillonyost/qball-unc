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
// ApplyElectricFieldPulseCmd.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef APPLYELECTRICFIELDPULSECMD_H
#define APPLYELECTRICFIELDPULSECMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"

class ApplyElectricFieldPulseCmd : public Cmd
{

  public:
  Sample *s;

  ApplyElectricFieldPulseCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "apply_electric_field_pulse"; }

  char const*help_msg(void) const
  {
    return
    "\n apply_electric_field_pulse e_direction e_strength\n\n"
    " syntax: apply_electric_field_pulse\n\n"
    "   The apply_electric_field_pulse apply an impulse electric field \n"
    " on the wavefunctions\n\n";
  }

  int action(int argc, char **argv);

  ApplyElectricFieldPulseCmd();
  ~ApplyElectricFieldPulseCmd();
};
#endif
