/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR VOLUME
/*
Calculate the volume of the simulation box.

\par Examples
The following input tells plumed to print the volume of the system
\verbatim
VOLUME LABEL=vol
PRINT ARG=vol
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC


class Volume : public Colvar {
  bool components;

public:
  Volume(const ActionOptions&);
// active methods:
  virtual void calculate();
/// Register all the keywords for this action
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Volume,"VOLUME")

Volume::Volume(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
components(false)
{
  std::vector<AtomNumber> atoms;
  parseFlag("COMPONENTS",components);
  checkRead();

  if(components){
// todo
  }
  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
}

void Volume::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); 
  keys.addFlag("COMPONENTS",false,"use xx, yy, zz, alpha, beta, gamma as the colvars rather than the box volume");
}


// calculator
void Volume::calculate(){
  if(components){
// todo
  };

  setBoxDerivatives(-1.0*Tensor::identity());
  setValue         (getBox().determinant());
}

}
}


