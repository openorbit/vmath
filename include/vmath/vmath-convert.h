/*
  Copyright 2006 Mattias Holm <mattias.holm(at)openorbit.org>

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


#ifndef MATH_CONVERT_H__
#define MATH_CONVERT_H__

#include <vmath/vmath-constants.h>
#include <vmath/vmath-integer.h>
#include <vmath/vmath-matvec.h>

#define DEG_TO_RAD(d) ((d) * M_PI/180.0f)
#define RAD_TO_DEG(r) ((r) * 180.0f/M_PI)

static inline double3
vf3_to_vd3(float3 v)
{
  return vd3_set(v.x, v.y, v.z);
}

static inline float3
vd3_to_vf3(double3 v)
{
  return vf3_set(v.x, v.y, v.z);
}

static inline float4
vd4_to_vf4(double4 v)
{
  return vf4_set(v.x, v.y, v.z, v.w);
}


static inline float3
v3i_to_v3f(int3 iv)
{
  float3 res = vf3_set((float)v3i_get(iv, 0),
                       (float)v3i_get(iv, 1),
                       (float)v3i_get(iv, 2));
  return res;
}

static inline double3
v3i_to_v3d(int3 iv)
{
  double3 res = vd3_set((double)v3i_get(iv, 0),
                        (double)v3i_get(iv, 1),
                        (double)v3i_get(iv, 2));
  return res;
}



static inline double3
v3l_to_v3d(long3 lv)
{
  double3 res = vd3_set((double)v3l_get(lv, 0),
                        (double)v3l_get(lv, 1),
                        (double)v3l_get(lv, 2));
  return res;
}

static inline float3
v3l_to_v3f(long3 lv)
{
  float3 res = vf3_set((float)v3l_get(lv, 0),
                       (float)v3l_get(lv, 1),
                       (float)v3l_get(lv, 2));
  return res;
}


static inline void
md4_to_mf4(const double4x4 d, float4x4 f)
{
  for (int i = 0 ; i < 4 ; i ++) {
    f[i] = vd4_to_vf4(d[i]);
  }
}



float3 equ2cart_f(float ra, float dec);
double3 equ2cart_d(double ra, double dec);
void
cart2geodetic_f(float3 p, float a, float e,
           float * restrict latitude,
           float * restrict longitude,
           float * restrict altitude);

float3
geodetic2cart_f(float a, float e,
                float latitude,
                float longitude,
                float altitude);

double3
geodetic2cart_d(double a, double e,
                double latitude,
                double longitude,
                double altitude);

#endif /* ! MATH_CONVERT_H__ */