/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2021 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

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
#include "tools/Torsion.h"

#include <string>
#include <cmath>

#include <iostream>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR TORSION
/*

*/
//+ENDPLUMEDOC

class MultiTorsion : public Colvar {
  unsigned maxn=3;
  unsigned ntorsions=0;
  std::vector<Value*> components;

public:
  explicit MultiTorsion(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(MultiTorsion,"MULTI_TORSION")

void MultiTorsion::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords( keys );
  keys.add("numbered","ATOMS","");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","NAME","");
  keys.add("compulsory","MAXN","3","");
  useCustomisableComponents(keys);
}

MultiTorsion::MultiTorsion(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  std::vector<AtomNumber> atoms;
  std::vector<std::string> names;
  maxn=3;
  parse("MAXN",maxn);
  for(int i=1;; ++i ) {
    std::vector<AtomNumber> t;
    parseAtomList("ATOMS", i, t );
    if( t.empty() ) break;
    if( t.size()!=4 ) {
      std::string ss; Tools::convert(i,ss);
      error("ATOMS" + ss + " keyword has the wrong number of atoms");
    }
    atoms.push_back(t[0]);
    atoms.push_back(t[1]);
    atoms.push_back(t[2]);
    atoms.push_back(t[3]);
    std::string name;
    parseNumbered("NAME",i,name);
    names.push_back(name);
    for(unsigned n=1; n<=maxn; n++) {
      std::string nstr;
      Tools::convert(n,nstr);
      addComponentWithDerivatives(name+"_sin_"+nstr); componentIsNotPeriodic(name+"_sin_"+nstr);
      components.push_back(getPntrToComponent(name+"_sin_"+nstr));
      addComponentWithDerivatives(name+"_cos_"+nstr); componentIsNotPeriodic(name+"_cos_"+nstr);
      components.push_back(getPntrToComponent(name+"_cos_"+nstr));
    }
  }

  requestAtoms(atoms);
  ntorsions=names.size();
}

// calculator
void MultiTorsion::calculate() {

  makeWhole();
  int k=0;
  for(unsigned i=0; i<ntorsions; i++) {
    Vector d0,d1,d2;
    d0=delta(getPosition(4*i+1),getPosition(4*i+0));
    d1=delta(getPosition(4*i+2),getPosition(4*i+1));
    d2=delta(getPosition(4*i+3),getPosition(4*i+2));
    double s,c;
    Vector ds_d0,ds_d1,ds_d2;
    Vector dc_d0,dc_d1,dc_d2;
    Torsion t;
    t.compute(d0,d1,d2,s,c,ds_d0,ds_d1,ds_d2,dc_d0,dc_d1,dc_d2);
    const double invR2=1.0/(c*c+s*s);
    const double sqr=std::sqrt(invR2);
    auto ss=s*sqr;
    auto cc=c*sqr;
// d sqrt(1/(s**2+c**2)) /ds =  -s/(s**2+c**2)**3/2= -s*sqr**3
    auto fs=-s*(sqr*sqr*sqr);
    auto fc=-c*(sqr*sqr*sqr);
    auto dss_d0=ds_d0*sqr+s*ds_d0*fs+s*dc_d0*fc;
    auto dss_d1=ds_d1*sqr+s*ds_d1*fs+s*dc_d1*fc;
    auto dss_d2=ds_d2*sqr+s*ds_d2*fs+s*dc_d2*fc;
    auto dcc_d0=dc_d0*sqr+c*dc_d0*fc+c*ds_d0*fs;
    auto dcc_d1=dc_d1*sqr+c*dc_d1*fc+c*ds_d1*fs;
    auto dcc_d2=dc_d2*sqr+c*dc_d2*fc+c*ds_d2*fs;
    // sin(nx)
    k=0;
    if(maxn>0) {
      components[k]->set(ss);
      setAtomsDerivatives(components[k],4*i+0,dss_d0);
      setAtomsDerivatives(components[k],4*i+1,-dss_d0+dss_d1);
      setAtomsDerivatives(components[k],4*i+2,-dss_d1+dss_d2);
      setAtomsDerivatives(components[k],4*i+3,-dss_d2);
      setBoxDerivativesNoPbc(components[k]);
      k++;
      components[k]->set(cc);
      setAtomsDerivatives(components[k],4*i+0,dcc_d0);
      setAtomsDerivatives(components[k],4*i+1,-dcc_d0+dcc_d1);
      setAtomsDerivatives(components[k],4*i+2,-dcc_d1+dcc_d2);
      setAtomsDerivatives(components[k],4*i+3,-dcc_d2);
      setBoxDerivativesNoPbc(components[k]);
      k++;
    }
    if(maxn>1) {
      auto ss2=2*ss*cc;
      auto dss2_d0=2*(ss*dcc_d0+dss_d0*cc);
      auto dss2_d1=2*(ss*dcc_d1+dss_d1*cc);
      auto dss2_d2=2*(ss*dcc_d2+dss_d2*cc);
      components[k]->set(ss2);
      setAtomsDerivatives(components[k],4*i+0,dss2_d0);
      setAtomsDerivatives(components[k],4*i+1,-dss2_d0+dss2_d1);
      setAtomsDerivatives(components[k],4*i+2,-dss2_d1+dss2_d2);
      setAtomsDerivatives(components[k],4*i+3,-dss2_d2);
      setBoxDerivativesNoPbc(components[k]);
      k++;
      auto cc2=1-2*ss*ss;
      auto dcc2_d0=-4*ss*dss_d0;
      auto dcc2_d1=-4*ss*dss_d1;
      auto dcc2_d2=-4*ss*dss_d2;
      components[k]->set(cc2);
      setAtomsDerivatives(components[k],4*i+0,dcc2_d0);
      setAtomsDerivatives(components[k],4*i+1,-dcc2_d0+dcc2_d1);
      setAtomsDerivatives(components[k],4*i+2,-dcc2_d1+dcc2_d2);
      setAtomsDerivatives(components[k],4*i+3,-dcc2_d2);
      setBoxDerivativesNoPbc(components[k]);
      k++;
    }
    if(maxn>2) {
      auto ss3=3*ss-4*ss*ss*ss;
      auto dss3_d0=(3-12*ss*ss)*dss_d0;
      auto dss3_d1=(3-12*ss*ss)*dss_d1;
      auto dss3_d2=(3-12*ss*ss)*dss_d2;
      components[k]->set(ss3);
      setAtomsDerivatives(components[k],4*i+0,dss3_d0);
      setAtomsDerivatives(components[k],4*i+1,-dss3_d0+dss3_d1);
      setAtomsDerivatives(components[k],4*i+2,-dss3_d1+dss3_d2);
      setAtomsDerivatives(components[k],4*i+3,-dss3_d2);
      setBoxDerivativesNoPbc(components[k]);
      k++;
      auto cc3=4*cc*cc*cc-3*cc;
      auto dcc3_d0=(12*cc*cc-3)*dcc_d0;
      auto dcc3_d1=(12*cc*cc-3)*dcc_d1;
      auto dcc3_d2=(12*cc*cc-3)*dcc_d2;
      components[k]->set(cc3);
      setAtomsDerivatives(components[k],4*i+0,dcc3_d0);
      setAtomsDerivatives(components[k],4*i+1,-dcc3_d0+dcc3_d1);
      setAtomsDerivatives(components[k],4*i+2,-dcc3_d1+dcc3_d2);
      setAtomsDerivatives(components[k],4*i+3,-dcc3_d2);
      setBoxDerivativesNoPbc(components[k]);
      k++;
    }
    if(maxn>3) plumed_error();
  }

}

}
}



