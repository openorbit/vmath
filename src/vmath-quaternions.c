/*
 Copyright 2006,2013 Mattias Holm <lorrden(at)openorbit.org>
 
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



#include <tgmath.h>

#include <vmath/vmath-types.h>
#include <vmath/vmath-constants.h>
#include <vmath/vmath-matvec.h>
#include <vmath/vmath-quaternions.h>
#include <assert.h>

float
qf_scalar(quatf_t q)
{
#if __has_feature(attribute_ext_vector_type)
  return q.w;
#else
  #error "not implemented"
#endif
}

double
qd_scalar(quatd_t q)
{
#if __has_feature(attribute_ext_vector_type)
  return q.w;
#else
  #error "not implemented"
#endif
}


float3
vf3_qf_rot(float3 v, quatf_t q)
{
  float3x3 m;
  qf_mf3_convert(m, q);

  float3 res = mf3_v_mul(m, v);
  return res;
}

double3
vd3_qd_rot(double3 v, quatd_t q)
{
  double3x3 m;
  qd_md3_convert(m, q);

  double3 res = md3_v_mul(m, v);
  return res;
}


float3
qf_vector(quatf_t q)
{
#if __has_feature(attribute_ext_vector_type)
  float3 r = {q.x, q.y, q.z};
  return r;
#else
#error "not implemented"
#endif
}

double3
qd_vector(quatd_t q)
{
#if __has_feature(attribute_ext_vector_type)
  double3 r = {q.x, q.y, q.z};
  return r;
#else
#error "not implemented"
#endif
}


void
qf_mf3_convert(float3x3 m, quatf_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  float n = qf_dot(q, q);
  float a = (n > 0.0f) ? 2.0f / n : 0.0f;

  float xa = q.x*a, ya = q.y*a, za = q.z*a;
  float xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  float yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  float wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vf3_set(1.0f-(yy+zz), xy-wz, xz+wy);
  m[1] = vf3_set(xy+wz, 1.0f-(xx+zz), yz-wx);
  m[2] = vf3_set(xz-wy, yz+wx, 1.0f-(xx+yy));
}


void
qd_md3_convert(double3x3 m, quatd_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  double n = qd_dot(q, q);
  double a = (n > 0.0) ? 2.0 / n : 0.0;

  double xa = q.x*a, ya = q.y*a, za = q.z*a;
  double xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  double yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  double wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vd3_set(1.0-(yy+zz), xy-wz, xz+wy);
  m[1] = vd3_set(xy+wz, 1.0-(xx+zz), yz-wx);
  m[2] = vd3_set(xz-wy, yz+wx, 1.0-(xx+yy));
}



void
qf_mf4_convert(float4x4 m, quatf_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  float n = qf_dot(q, q);
  float a = (n > 0.0f) ? 2.0f / n : 0.0f;

  float xa = q.x*a, ya = q.y*a, za = q.z*a;
  float xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  float yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  float wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vf4_set(1.0f-(yy+zz), xy-wz, xz+wy, 0.0f);
  m[1] = vf4_set(xy+wz, 1.0f-(xx+zz), yz-wx, 0.0f);
  m[2] = vf4_set(xz-wy, yz+wx, 1.0f-(xx+yy), 0.0f);
  m[3] = vf4_set(0.0f, 0.0f, 0.0f, 1.0f);
}


void
qd_md4_convert(double4x4 m, quatd_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  double n = qd_dot(q, q);
  double a = (n > 0.0) ? 2.0 / n : 0.0;

  double xa = q.x*a, ya = q.y*a, za = q.z*a;
  double xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  double yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  double wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vd4_set(1.0-(yy+zz), xy-wz, xz+wy, 0.0);
  m[1] = vd4_set(xy+wz, 1.0-(xx+zz), yz-wx, 0.0);
  m[2] = vd4_set(xz-wy, yz+wx, 1.0-(xx+yy), 0.0);
  m[3] = vd4_set(0.0, 0.0, 0.0, 1.0);
}


void
qf_mf4_convert_inv(float4x4 m, quatf_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  float n = qf_dot(q, q);
  float a = (n > 0.0f) ? 2.0f / n : 0.0f;

  float xa = q.x*a, ya = q.y*a, za = q.z*a;
  float xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  float yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  float wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vf4_set(1.0f-(yy+zz), xy+wz, xz-wy, 0.0f);
  m[1] = vf4_set(xy-wz, 1.0f-(xx+zz), yz+wx, 0.0f);
  m[2] = vf4_set(xz+wy, yz-wx, 1.0f-(xx+yy), 0.0f);
  m[3] = vf4_set(0.0f, 0.0f, 0.0f, 1.0f);
}

void
qd_md4_convert_inv(double4x4 m, quatd_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  double n = qd_dot(q, q);
  double a = (n > 0.0f) ? 2.0f / n : 0.0f;

  double xa = q.x*a, ya = q.y*a, za = q.z*a;
  double xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  double yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  double wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vd4_set(1.0-(yy+zz), xy+wz, xz-wy, 0.0);
  m[1] = vd4_set(xy-wz, 1.0-(xx+zz), yz+wx, 0.0);
  m[2] = vd4_set(xz+wy, yz-wx, 1.0-(xx+yy), 0.0);
  m[3] = vd4_set(0.0, 0.0, 0.0, 1.0);
}


void
qf_mf3_convert_inv(float3x3 m, quatf_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  float n = qf_dot(q, q);
  float a = (n > 0.0f) ? 2.0f / n : 0.0f;

  float xa = q.x*a, ya = q.y*a, za = q.z*a;
  float xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  float yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  float wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vf3_set(1.0f-(yy+zz), xy+wz, xz-wy);
  m[1] = vf3_set(xy-wz, 1.0f-(xx+zz), yz+wx);
  m[2] = vf3_set(xz+wy, yz-wx, 1.0f-(xx+yy));
}

void
qd_md3_convert_inv(double3x3 m, quatd_t q)
{
#if ! __has_feature(attribute_ext_vector_type)
#error "clang extended vector attributes required"
#endif

  double n = qd_dot(q, q);
  double a = (n > 0.0) ? 2.0 / n : 0.0;

  double xa = q.x*a, ya = q.y*a, za = q.z*a;
  double xx = q.x*xa, xy = q.x*ya, xz = q.x*za;
  double yy = q.y*ya, yz = q.y*za, zz = q.z*za;
  double wx = q.w*xa, wy = q.w*ya, wz = q.w*za;

  m[0] = vd3_set(1.0-(yy+zz), xy+wz, xz-wy);
  m[1] = vd3_set(xy-wz, 1.0-(xx+zz), yz+wx);
  m[2] = vd3_set(xz+wy, yz-wx, 1.0-(xx+yy));
}


quatf_t
qf_slerp(quatf_t q0, quatf_t q1, float t)
{
  // See http://en.wikipedia.org/wiki/Slerp
  if (t >= 1.0) return q1;
  if (t <= 0.0) return q0;

  float qdot = qf_dot(q0, q1);

  quatf_t q1prim;
  if (qdot < 0.0) {
    q1prim = -q1;
    qdot = -qdot;
  } else {
    q1prim = q1;
  }

  if (qdot < -1.0)  qdot = -1.0;
  if (qdot > 1.0)  qdot = 1.0;
  assert(qdot >= -1.0);
  assert(qdot <= 1.0);

  float qang = acos(qdot);
  float s0 = sin((1.0-t)*qang) / sin(qang);
  float s1 = sin(t*qang) / sin(qang);

  if (qang) { // If slerping between the same points, qang will be 0, so s0 will be NaN or Inf.
    quatf_t res = s0 * q0 + s1 * q1prim;
    return res;
  }

  return q0;
}

quatd_t
qd_slerp(quatd_t q0, quatd_t q1, double t)
{
  // See http://en.wikipedia.org/wiki/Slerp
  if (t >= 1.0) return q1;
  if (t <= 0.0) return q0;

  double qdot = qd_dot(q0, q1);

  quatd_t q1prim;
  if (qdot < 0.0) {
    q1prim = -q1;
    qdot = -qdot;
  } else {
    q1prim = q1;
  }

  if (qdot < -1.0)  qdot = -1.0;
  if (qdot > 1.0)  qdot = 1.0;
  assert(qdot >= -1.0);
  assert(qdot <= 1.0);

  double qang = acos(qdot);
  double s0 = sin((1.0-t)*qang) / sin(qang);
  double s1 = sin(t*qang) / sin(qang);

  if (qang) { // If slerping between the same points, qang will be 0, so s0 will be NaN or Inf.
    quatd_t res = s0 * q0 + s1 * q1prim;
    return res;
  }

  return q0;
}



quatf_t
mf3_qf_convert(float3x3 m)
{
  quatf_t q;
  float tr, s;

  tr = m[0][0] + m[1][1] + m[2][2];

  if (tr >= 0.0f) {
    s = sqrtf(tr+1.0f);
    q.w = s*0.5f;
    s = 0.5f / s;
    q.x = (m[2][1] - m[1][2]) * s;
    q.y = (m[0][2] - m[2][0]) * s;
    q.z = (m[1][0] - m[0][1]) * s;
  } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
    s =  2.0f * sqrtf(m[0][0] - m[1][1] - m[2][2] + 1.0f);
    q.x = 0.25f * s;
    q.y = (m[0][1] + m[1][0] ) / s;
    q.z = (m[0][2] + m[2][0] ) / s;
    q.w = (m[2][1] - m[1][2] ) / s;
  } else if (m[1][1] > m[2][2]) {
    s = 2.0f * sqrtf(m[1][1] - m[0][0] - m[2][2] + 1.0f);
    q.x = (m[0][1] + m[1][0] ) / s;
    q.y = 0.25f * s;
    q.z = (m[1][2] + m[2][1] ) / s;
    q.w = (m[0][2] - m[2][0] ) / s;
  } else {
    s = 2.0f * sqrtf(m[2][2] - m[0][0] - m[1][1] + 1.0f);
    q.x = (m[0][2] + m[2][0]) / s;
    q.y = (m[1][2] + m[2][1]) / s;
    q.z = 0.25f * s;
    q.w = (m[1][0] - m[0][1]) / s;
  }

  return q;
}

quatd_t
md3_qd_convert(double3x3 m)
{
  quatd_t q;
  double tr, s;

  tr = m[0][0] + m[1][1] + m[2][2];

  if (tr >= 0.0f) {
    s = sqrtf(tr+1.0);
    q.w = s*0.5;
    s = 0.5 / s;
    q.x = (m[2][1] - m[1][2]) * s;
    q.y = (m[0][2] - m[2][0]) * s;
    q.z = (m[1][0] - m[0][1]) * s;
  } else if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
    s =  2.0 * sqrt(m[0][0] - m[1][1] - m[2][2] + 1.0);
    q.x = 0.25 * s;
    q.y = (m[0][1] + m[1][0] ) / s;
    q.z = (m[0][2] + m[2][0] ) / s;
    q.w = (m[2][1] - m[1][2] ) / s;
  } else if (m[1][1] > m[2][2]) {
    s = 2.0 * sqrt(m[1][1] - m[0][0] - m[2][2] + 1.0);
    q.x = (m[0][1] + m[1][0] ) / s;
    q.y = 0.25f * s;
    q.z = (m[1][2] + m[2][1] ) / s;
    q.w = (m[0][2] - m[2][0] ) / s;
  } else {
    s = 2.0 * sqrt(m[2][2] - m[0][0] - m[1][1] + 1.0);
    q.x = (m[0][2] + m[2][0]) / s;
    q.y = (m[1][2] + m[2][1]) / s;
    q.z = 0.25f * s;
    q.w = (m[1][0] - m[0][1]) / s;
  }

  return q;
}



quatf_t
mf4_qf_convert(float4x4 m)
{
  quatf_t q;

  float tr, s;
  tr = m[0][0] + m[1][1] + m[2][2];
  if (tr >= 0.0f) {
    s = sqrt(tr+m[3][3]);
    q.w = s*0.5;
    s = 0.5 / s;
    q.x = (m[2][1] - m[1][2]) * s;
    q.y = (m[0][2] - m[2][0]) * s;
    q.z = (m[1][0] - m[0][1]) * s;
  } else {
    int h = 0;
    if (m[1][1] > m[0][0]) h = 1;
    if (m[2][2] > m[h][h]) h = 2;
    switch (h) {
#define CASE_MACRO(i,j,k,I,J,K)                             \
  case I:                                                   \
    s = sqrt( (m[I][I] - (m[J][J]+m[K][K])) + m[3][3] );    \
    q[i] = s*0.5;                                           \
    q[j] = (m[I][J] + m[J][I]) * s;                         \
    q[k] = (m[K][I] + m[I][K]) * s;                         \
    q.w = (m[K][J] - m[J][K]) * s;                           \
    break

        CASE_MACRO(0, 1, 2, 0, 1, 2);
        CASE_MACRO(1, 2, 0, 1, 2, 0);
        CASE_MACRO(2, 0, 1, 2, 0, 1);
#undef CASE_MACRO
      default:
        assert(0);
    }
  }

  // QUERY: Is the last ref to z correct?
  if (m[3][3] != 0.0f) {
    s = 1.0f / sqrt(m[3][3]);
    q.x *= s; q.y *= s; q.z *= s; //QZ(q) *= s;
  }
  return q;
}

quatd_t
md4_qd_convert(double4x4 m)
{
  quatd_t q;

  double tr, s;
  tr = m[0][0] + m[1][1] + m[2][2];
  if (tr >= 0.0f) {
    s = sqrt(tr+m[3][3]);
    q.w = s*0.5;
    s = 0.5 / s;
    q.x = (m[2][1] - m[1][2]) * s;
    q.y = (m[0][2] - m[2][0]) * s;
    q.z = (m[1][0] - m[0][1]) * s;
  } else {
    int h = 0;
    if (m[1][1] > m[0][0]) h = 1;
    if (m[2][2] > m[h][h]) h = 2;
    switch (h) {
#define CASE_MACRO(i,j,k,I,J,K)                           \
  case I:                                                 \
    s = sqrt( (m[I][I] - (m[J][J]+m[K][K])) + m[3][3] );  \
    q[i] = s*0.5;                                         \
    q[j] = (m[I][J] + m[J][I]) * s;                       \
    q[k] = (m[K][I] + m[I][K]) * s;                       \
    q.w = (m[K][J] - m[J][K]) * s;                        \
    break

    CASE_MACRO(0, 1, 2, 0, 1, 2);
    CASE_MACRO(1, 2, 0, 1, 2, 0);
    CASE_MACRO(2, 0, 1, 2, 0, 1);
#undef CASE_MACRO
    default:
      assert(0);
    }
  }

  // QUERY: Is the last ref to z correct?
  if (m[3][3] != 0.0f) {
    s = 1.0f / sqrt(m[3][3]);
    q.x *= s; q.y *= s; q.z *= s; //QZ(q) *= s;
  }
  return q;
}




quatf_t
qf_add(const quatf_t a, const quatf_t b)
{
  return vf4_add(a, b);
}

quatd_t
qd_add(const quatd_t a, const quatd_t b)
{
  return vd4_add(a, b);
}


quatf_t
qf_mul(const quatf_t a, const quatf_t b)
{
#if __has_feature(attribute_ext_vector_type)
  quatf_t r;
  r.x = a.x*b.w + a.w*b.x	+ a.y*b.z - a.z*b.y;
  r.y = a.y*b.w + a.w*b.y	+ a.z*b.x - a.x*b.z;
  r.z = a.z*b.w + a.w*b.z	+ a.x*b.y - a.y*b.x;
  r.w = a.w*b.w - a.x*b.x	- a.y*b.y - a.z*b.z;
  return r;
#else
  #error "not implemented"
#endif
}

quatd_t
qd_mul(const quatd_t a, const quatd_t b)
{
#if __has_feature(attribute_ext_vector_type)
  quatd_t r;
  r.x = a.x*b.w + a.w*b.x	+ a.y*b.z - a.z*b.y;
  r.y = a.y*b.w + a.w*b.y	+ a.z*b.x - a.x*b.z;
  r.z = a.z*b.w + a.w*b.z	+ a.x*b.y - a.y*b.x;
  r.w = a.w*b.w - a.x*b.x	- a.y*b.y - a.z*b.z;
  return r;
#else
  #error "not implemented"
#endif
}


quatf_t
qf_s_div(quatf_t q, float d)
{
#if __has_feature(attribute_ext_vector_type)
  quatf_t r;
  r.x = q.x / d;
  r.y = q.y / d;
  r.z = q.z / d;
  r.w = q.w / d;
  return r;
#else
  #error "not implemented"
#endif
}

quatd_t
qd_s_div(quatd_t q, double d)
{
#if __has_feature(attribute_ext_vector_type)
  quatd_t r;
  r.x = q.x / d;
  r.y = q.y / d;
  r.z = q.z / d;
  r.w = q.w / d;
  return r;
#else
#error "not implemented"
#endif
}


quatf_t
qf_s_mul(quatf_t q, float d)
{
#if __has_feature(attribute_ext_vector_type)
  quatf_t r;
  r = q * d;
  return r;
#else
#error "not implemented"
#endif
}

quatd_t
qd_s_mul(quatd_t q, double d)
{
#if __has_feature(attribute_ext_vector_type)
  quatd_t r;
  r = q * d;
  return r;
#else
#error "not implemented"
#endif
}




float
qf_dot(quatf_t a, quatf_t b)
{
  return vf4_dot(a, b);
}

double
qd_dot(quatd_t a, quatd_t b)
{
  return vd4_dot(a, b);
}

float3
qf_cross(const quatf_t a, const quatf_t b)
{
  return vf3_cross(qf_vector(a), qf_vector(b));
}

double3
qd_cross(const quatd_t a, const quatd_t b)
{
  return vd3_cross(qd_vector(a), qd_vector(b));
}


float
qf_abs(quatf_t q)
{
  return vf4_abs(q);
}

double
qd_abs(quatd_t q)
{
  return vd4_abs(q);
}


quatf_t
qf_conj(const quatf_t q)
{
#if __has_feature(attribute_ext_vector_type)
  quatf_t qp = {-q.x, -q.y, -q.z, q.w};
  return qp;
#else
  #error "not implemented"
#endif
}

quatd_t
qd_conj(const quatd_t q)
{
#if __has_feature(attribute_ext_vector_type)
  quatd_t qp = {-q.x, -q.y, -q.z, q.w};
  return qp;
#else
#error "not implemented"
#endif
}

quatf_t
qf_repr(const quatf_t q)
{
  quatf_t res;
  quatf_t qp;
  qp = qf_conj(q);
  float d = qf_dot(q, q);
  res = qf_s_div(qp, d);
  return res;
}

quatd_t
qd_repr(const quatd_t q)
{
  quatd_t res;
  quatd_t qp;
  qp = qd_conj(q);
  double d = qd_dot(q, q);
  res = qd_s_div(qp, d);
  return res;
}


quatf_t
qf_div(quatf_t a, quatf_t b)
{
  quatf_t res;
  quatf_t br;

  br = qf_repr(b);

  res = qf_mul(a, br);
  return res;
}


quatd_t
qd_div(quatd_t a, quatd_t b)
{
  quatd_t res;
  quatd_t br;

  br = qd_repr(b);

  res = qd_mul(a, br);
  return res;
}



quatf_t
qf_rotv(float3 axis, float alpha)
{
  quatf_t q;
  float Omega = alpha * 0.5;
  float sin_Omega = sin(Omega);
#if __has_feature(attribute_ext_vector_type)
  q.x = axis.x * sin_Omega;
  q.y = axis.y * sin_Omega;
  q.z = axis.z * sin_Omega;
  q.w = cos(Omega);
#else
#error "not implemented"
#endif
  return q;
}

quatd_t
qd_rotv(double3 axis, double alpha)
{
  quatd_t q;
  double Omega = alpha * 0.5;
  double sin_Omega = sin(Omega);
#if __has_feature(attribute_ext_vector_type)
  q.x = axis.x * sin_Omega;
  q.y = axis.y * sin_Omega;
  q.z = axis.z * sin_Omega;
  q.w = cos(Omega);
#else
#error "not implemented"
#endif
  return q;
}


quatf_t
qf_rot(float x, float y, float z, float alpha)
{
  float3 axis = {x, y, z};
  return qf_rotv(axis, alpha);
}

quatd_t
qd_rot(double x, double y, double z, double alpha)
{
  double3 axis = {x, y, z};
  return qd_rotv(axis, alpha);
}

quatf_t
qf_normalise(quatf_t q)
{
  return vf4_normalise(q);
}

quatd_t
qd_normalise(quatd_t q)
{
  return vd4_normalise(q);
}
