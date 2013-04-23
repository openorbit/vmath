/*
 Copyright 2013 Mattias Holm <lorrden(at)openorbit.org>

 This file is part of Open Orbit. Open Orbit is free software: you can
 redistribute it and/or modify it under the terms of the GNU General Public
 License as published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 You should have received a copy of the GNU General Public License
 along with Open Orbit.  If not, see <http://www.gnu.org/licenses/>.

 Some files of Open Orbit have relaxed licensing conditions. This file is
 licenced under the 2-clause BSD licence.

 Redistribution and use of this file in source and binary forms, with or
 without modification, are permitted provided that the following conditions are
 met:

 - Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 - Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <vmath/vmath.h>
#include <math.h>
#include <assert.h>

float3
vf3_cube_sphere_map(float3 p)
{
  float3 tmp;
  tmp.x = sqrtf(1.0-p.y*p.y*0.5-p.z*p.z*0.5+p.y*p.y*p.z*p.z/3.0);
  tmp.y = sqrtf(1.0-p.z*p.z*0.5-p.x*p.x*0.5+p.z*p.z*p.x*p.x/3.0);
  tmp.z = sqrtf(1.0-p.x*p.x*0.5-p.y*p.y*0.5+p.x*p.x*p.y*p.y/3.0);

  float3 res = vf3_mul(p, tmp);
  return res;
}

// The function is not hardened enough at the moment. Doubles are used, as
// rounding errors will otherwise result in negative 0 to sqrt.
float3
vf3_sphere_cube_map(float3 p)
{
  float3 pabs = vf3_set(fabsf(p.x), fabsf(p.y), fabsf(p.z));
  float3 cp;

  if (pabs.x >= pabs.y && pabs.x >= pabs.z) {
    cp.x = copysign(1.0, p.x);

    double y2 = 2.0*p.y*p.y;
    double z2 = 2.0*p.z*p.z;
    double y2z2_3 = -y2+z2-3.0;
    double inner = -sqrt(y2z2_3*y2z2_3 - 12.0*y2);

    if (p.y == 0.0) cp.y = 0.0;
    else cp.y = copysign(sqrt(inner + y2 - z2 + 3.0)/M_SQRT2, p.y);

    if (p.z == 0.0) cp.z = 0.0;
    else cp.z = copysign(sqrt(inner - y2 + z2 + 3.0)/M_SQRT2, p.z);
  } else if (pabs.y >= pabs.x && pabs.y >= pabs.z) {
    double z2 = 2.0*p.z*p.z;
    double x2 = 2.0*p.x*p.x;
    double z2x2_3 = -z2+x2-3.0;
    double inner = -sqrt(z2x2_3*z2x2_3-12.0*z2);

    if (p.x == 0.0) cp.x = 0.0;
    else cp.x = copysign(sqrt(inner + x2 - z2 + 3.0)/M_SQRT2, p.x);

    cp.y = copysign(1.0, p.y);

    if (p.z == 0.0) cp.z = 0.0;
    else cp.z = copysign(sqrt(inner - x2 + z2 + 3.0)/M_SQRT2, p.z);
  } else if (pabs.z >= pabs.x && pabs.z >= pabs.y) {
    double x2 = 2.0*p.x*p.x;
    double y2 = 2.0*p.y*p.y;
    double x2y2_3 = -x2+y2-3.0;
    double inner = -sqrt(x2y2_3*x2y2_3-12.0*x2);

    if (p.x == 0.0) cp.x = 0.0;
    else cp.x = copysign(sqrt(inner + x2 - y2 + 3.0)/M_SQRT2, p.x);

    if (p.y == 0.0) cp.y = 0.0;
    else cp.y = copysign(sqrt(inner - x2 + y2 + 3.0)/M_SQRT2, p.y);

    cp.z = copysign(1.0, p.z);
  } else {
    assert(0 && "impossible");
  }

  return cp;
}


double3
vd3_cube_sphere_map(double3 p)
{
  double3 tmp;
  tmp.x = sqrtf(1.0-p.y*p.y*0.5-p.z*p.z*0.5+p.y*p.y*p.z*p.z/3.0);
  tmp.y = sqrtf(1.0-p.z*p.z*0.5-p.x*p.x*0.5+p.z*p.z*p.x*p.x/3.0);
  tmp.z = sqrtf(1.0-p.x*p.x*0.5-p.y*p.y*0.5+p.x*p.x*p.y*p.y/3.0);

  double3 res = vd3_mul(p, tmp);
  return res;
}

// The function is not hardened enough at the moment. Doubles are used, as
// rounding errors will otherwise result in negative 0 to sqrt.
double3
vd3_sphere_cube_map(double3 p)
{
  double3 pabs = vd3_set(fabs(p.x), fabs(p.y), fabs(p.z));
  double3 cp;

  if (pabs.x >= pabs.y && pabs.x >= pabs.z) {
    cp.x = copysign(1.0, p.x);

    double y2 = 2.0*p.y*p.y;
    double z2 = 2.0*p.z*p.z;
    double y2z2_3 = -y2+z2-3.0;
    double inner = -sqrt(y2z2_3*y2z2_3 - 12.0*y2);

    if (p.y == 0.0) cp.y = 0.0;
    else cp.y = copysign(sqrt(inner + y2 - z2 + 3.0)/M_SQRT2, p.y);

    if (p.z == 0.0) cp.z = 0.0;
    else cp.z = copysign(sqrt(inner - y2 + z2 + 3.0)/M_SQRT2, p.z);
  } else if (pabs.y >= pabs.x && pabs.y >= pabs.z) {
    double z2 = 2.0*p.z*p.z;
    double x2 = 2.0*p.x*p.x;
    double z2x2_3 = -z2+x2-3.0;
    double inner = -sqrt(z2x2_3*z2x2_3-12.0*z2);

    if (p.x == 0.0) cp.x = 0.0;
    else cp.x = copysign(sqrt(inner + x2 - z2 + 3.0)/M_SQRT2, p.x);

    cp.y = copysign(1.0, p.y);

    if (p.z == 0.0) cp.z = 0.0;
    else cp.z = copysign(sqrt(inner - x2 + z2 + 3.0)/M_SQRT2, p.z);
  } else if (pabs.z >= pabs.x && pabs.z >= pabs.y) {
    double x2 = 2.0*p.x*p.x;
    double y2 = 2.0*p.y*p.y;
    double x2y2_3 = -x2+y2-3.0;
    double inner = -sqrt(x2y2_3*x2y2_3-12.0*x2);

    if (p.x == 0.0) cp.x = 0.0;
    else cp.x = copysign(sqrt(inner + x2 - y2 + 3.0)/M_SQRT2, p.x);

    if (p.y == 0.0) cp.y = 0.0;
    else cp.y = copysign(sqrt(inner - x2 + y2 + 3.0)/M_SQRT2, p.y);

    cp.z = copysign(1.0, p.z);
  } else {
    assert(0 && "impossible");
  }
  
  return cp;
}

