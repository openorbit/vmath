/*
 Copyright 2009 Mattias Holm <mattias.holm(at)openorbit.org>

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
#include <vmath/lwcoord.h>
#include <vmath/vmath-convert.h>

#include <stdio.h>

#ifndef __has_feature
#define __has_feature(x) 0  // Compatibility with non-clang compilers.
#endif


#if __has_feature(attribute_ext_vector_type)
#define OFFS_X(c) ((c)->offs.x)
#define OFFS_Y(c) ((c)->offs.y)
#define OFFS_Z(c) ((c)->offs.z)

#define SEG_X(c) ((c)->seg.x)
#define SEG_Y(c) ((c)->seg.y)
#define SEG_Z(c) ((c)->seg.z)

#define OFFS_X64(c) ((c)->offs.x)
#define OFFS_Y64(c) ((c)->offs.y)
#define OFFS_Z64(c) ((c)->offs.z)

#define SEG_X64(c) ((c)->seg.x)
#define SEG_Y64(c) ((c)->seg.y)
#define SEG_Z64(c) ((c)->seg.z)

#else

#define OFFS_X(c) (((float*)&(c)->offs)[0])
#define OFFS_Y(c) (((float*)&(c)->offs)[1])
#define OFFS_Z(c) (((float*)&(c)->offs)[2])

#define SEG_X(c) (((int32_t*)&(c)->seg)[0])
#define SEG_Y(c) (((int32_t*)&(c)->seg)[1])
#define SEG_Z(c) (((int32_t*)&(c)->seg)[2])


#define OFFS_X64(c) (((double*)&(c)->offs)[0])
#define OFFS_Y64(c) (((double*)&(c)->offs)[1])
#define OFFS_Z64(c) (((double*)&(c)->offs)[2])

#define SEG_X64(c) (((int64_t*)&(c)->seg)[0])
#define SEG_Y64(c) (((int64_t*)&(c)->seg)[1])
#define SEG_Z64(c) (((int64_t*)&(c)->seg)[2])
#endif
void
lwc_set(lwcoord_t *coord, double x, double y, double z)
{
  coord->offs = vd3_set(x, y, z);
  coord->seg = vl3_set(0, 0, 0);

  lwc_normalise(coord);
}

void
lwc_setv(lwcoord_t *coord, double3 v)
{
  coord->offs = v;
  coord->seg = vl3_set(0, 0, 0);
  lwc_normalise(coord);
}


void
lwc_normalise(lwcoord_t *coord)
{
  if (fabsf(coord->offs.x) >= OO_LW_SEGMENT_LEN) {
    coord->seg.x += (int64_t) (coord->offs.x / OO_LW_SEGMENT_LEN);
    coord->offs.x = fmodf(coord->offs.x, OO_LW_SEGMENT_LEN);
  }
  if (fabsf(coord->offs.y) >= OO_LW_SEGMENT_LEN) {
    coord->seg.y += (int64_t) (coord->offs.y / OO_LW_SEGMENT_LEN);
    coord->offs.y = fmodf(coord->offs.y, OO_LW_SEGMENT_LEN);
  }
  if (fabsf(coord->offs.z) >= OO_LW_SEGMENT_LEN) {
    coord->seg.z += (int64_t) (coord->offs.z / OO_LW_SEGMENT_LEN);
    coord->offs.z = fmodf(coord->offs.z, OO_LW_SEGMENT_LEN);
  }
}

void
lwc_mul(lwcoord_t *lwc, float b)
{
  double3 p = lwc_global(lwc);
  p = p * b;
  lwc_set(lwc, p.x, p.y, p.z);
}

void
lwc_div(lwcoord_t *lwc, float b)
{
  double3 p = lwc_global(lwc);
  p = p / b;
  lwc_set(lwc, p.x, p.y, p.z);
}

void
lwc_dump(const lwcoord_t *lwc)
{
  fprintf(stderr, "lwc: [%lld %lld %lld]/[%f %f %f]\n",
          lwc->seg.x, lwc->seg.y, lwc->seg.z,
          lwc->offs.x, lwc->offs.y, lwc->offs.z);
}


void
lwc_translate3fv(lwcoord_t *coord, float3 offs)
{
  double3 offsd = vf3_to_vd3(offs);
  coord->offs = vd3_add(coord->offs, offsd);
  lwc_normalise(coord);
}

void
lwc_translate3dv(lwcoord_t *coord, double3 offs)
{
  double3 v = vd3_add(coord->offs, offs);
  coord->offs = v;
  lwc_normalise(coord);
}


void
lwc_translate3f(lwcoord_t *coord, float dx, float dy, float dz)
{
  coord->offs = vd3_add(coord->offs, vd3_set(dx, dy, dz));
  lwc_normalise(coord);
}

double3
lwc_globald(const lwcoord_t *coord)
{
  double3 p = coord->offs;
  double3 seg = v3l_to_v3d(coord->seg);

  return vd3_add(p, vd3_s_mul(seg, OO_LW_SEGMENT_LEN));
}


double3
lwc_global(const lwcoord_t *coord)
{
  double3 p = coord->offs;
  double3 seg = v3l_to_v3d(coord->seg);

  return vd3_add(p, vd3_s_mul(seg, OO_LW_SEGMENT_LEN));
}

double3
lwc_relvec(const lwcoord_t *coord, long3 seg)
{
  double3 r = coord->offs;
  long3 segdiff = coord->seg - seg;
  double3 segdiffr = vd3_set((double)segdiff.x,
                             (double)segdiff.y,
                             (double)segdiff.z);
  r = vd3_add(r, vd3_s_mul(segdiffr, OO_LW_SEGMENT_LEN));
  return r;
}

double3
lwc_relvec_d3(const lwcoord_t *coord, double3 center)
{
  double3 cordp = lwc_global(coord);

  return vd3_sub(cordp, center);
}

double3
lwc_dist(const lwcoord_t *a, const lwcoord_t * b)
{
  double3 diff = vd3_sub(a->offs, b->offs);
  long3 segdiff = a->seg - b->seg;
  double3 segdiffr = vd3_set((double)segdiff.x,
                             (double)segdiff.y,
                             (double)segdiff.z);

  return vd3_add(diff, vd3_s_mul(segdiffr, OO_LW_SEGMENT_LEN));
}


int
lwc_octant(const lwcoord_t *a, const lwcoord_t * b)
{
  double3 rel = lwc_dist(b, a);

  int octant = 0;
  octant = vd3_octant(vd3_set(0, 0, 0), rel);
  return octant;
}
