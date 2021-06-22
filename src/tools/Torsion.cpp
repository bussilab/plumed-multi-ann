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
#include "Torsion.h"
#include "Tensor.h"

#include <cmath>
#include <iostream>

namespace PLMD {

double Torsion::compute(const Vector& v1,const Vector& v2,const Vector& v3)const {
  const Vector nv2(v2*(1.0/v2.modulo()));
  const Vector a(crossProduct(nv2,v1));
  const Vector b(crossProduct(v3,nv2));
  const double cosangle=dotProduct(a,b);
  const double sinangle=dotProduct(crossProduct(a,b),nv2);
  return std::atan2(-sinangle,cosangle);
}

static void  compute_local(const Vector& v1,const Vector& v2,const Vector& v3,
                           double & sinangle,double &cosangle,
                           Vector & ds_dv1, Vector & ds_dv2, Vector & ds_dv3,
                           Vector & dc_dv1, Vector & dc_dv2, Vector & dc_dv3) {
  const double modv2(1./v2.modulo());
  const Vector nv2(v2*modv2);
  const Tensor dnv2_v2((Tensor::identity()-extProduct(nv2,nv2))*modv2);

  const Vector a(crossProduct(v2,v1));
  const Tensor da_dv2(dcrossDv1(v2,v1));
  const Tensor da_dv1(dcrossDv2(v2,v1));
  const Vector b(crossProduct(v3,v2));
  const Tensor db_dv3(dcrossDv1(v3,v2));
  const Tensor db_dv2(dcrossDv2(v3,v2));

  cosangle=dotProduct(a,b);
  dc_dv1=matmul(b,da_dv1);
  dc_dv2=matmul(b,da_dv2) + matmul(a,db_dv2);
  dc_dv3=matmul(a,db_dv3);

  const Vector cab(crossProduct(a,b));
  const Tensor dcab_dv1(matmul(dcrossDv1(a,b),da_dv1));
  const Tensor dcab_dv2(matmul(dcrossDv1(a,b),da_dv2) + matmul(dcrossDv2(a,b),db_dv2));
  const Tensor dcab_dv3(matmul(dcrossDv2(a,b),db_dv3));

  sinangle=-dotProduct(cab,nv2);
  ds_dv1=-matmul(nv2,dcab_dv1);
  ds_dv2=-matmul(nv2,dcab_dv2)-matmul(cab,dnv2_v2);
  ds_dv3=-matmul(nv2,dcab_dv3);
}

double Torsion::compute(const Vector& v1,const Vector& v2,const Vector& v3,Vector& d1,Vector& d2,Vector& d3)const {


  double sinangle,cosangle;
  Vector ds_dv1,ds_dv2,ds_dv3;
  Vector dc_dv1,dc_dv2,dc_dv3;
  compute_local(v1,v2,v3,sinangle,cosangle,ds_dv1,ds_dv2,ds_dv3,dc_dv1,dc_dv2,dc_dv3);
  const double torsion=std::atan2(sinangle,cosangle);
// this is required since v1 and v3 are not normalized:
  const double invR2=1.0/(cosangle*cosangle+sinangle*sinangle);

  d1= ( ds_dv1*cosangle - sinangle * dc_dv1 ) *invR2;
  d2= ( ds_dv2*cosangle - sinangle * dc_dv2 ) *invR2;
  d3= ( ds_dv3*cosangle - sinangle * dc_dv3 ) *invR2;

  return torsion;
}

void Torsion::compute(const Vector& v1,const Vector& v2,const Vector& v3,double&sinangle,double&cosangle,Vector& ds_d1,Vector& ds_d2,Vector& ds_d3,Vector& dc_d1,Vector& dc_d2,Vector& dc_d3)const {
  compute_local(v1,v2,v3,sinangle,cosangle,ds_d1,ds_d2,ds_d3,dc_d1,dc_d2,dc_d3);
}


}



