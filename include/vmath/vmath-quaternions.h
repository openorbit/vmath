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


#ifndef __QUATERNIONS_H__
#define __QUATERNIONS_H__
#ifdef __cplusplus
extern "C" {
#endif 


#include <vmath/vmath-types.h>
#include <vmath/vmath-constants.h>
#include <vmath/vmath-matvec.h>
    
float qf_scalar(const quatf_t q)
  __attribute__ ((__pure__));
double qd_scalar(const quatd_t q)
  __attribute__ ((__pure__));


float3 qf_vector(const quatf_t q);
double3 qd_vector(const quatd_t q);

float3 vf3_qf_rot(float3 v, quatf_t q);
double3 vd3_qd_rot(double3 v, quatd_t q);

#define QF_IDENT vf4_set(0.0, 0.0, 0.0, 1.0)
#define QD_IDENT vd4_set(0.0, 0.0, 0.0, 1.0)

/*!
 * \brief Converts a quaternion to a rotation matrix.
 * 
 * The q_m_convert function converts quaternion q into the rotational matrix m.
 * The rotational matrix is in a right handed coordinate system. And thus
 * compatible with toolkits such as OpenGL.
 * 
 * \param m A reference to the resulting matrix.
 * \param q The quaternion.
 */

void qf_mf3_convert(float3x3 m, quatf_t q);
void qd_md3_convert(double3x3 m, quatd_t q);
void qf_mf3_convert_inv(float3x3 m, quatf_t q);
void qd_md3_convert_inv(double3x3 m, quatd_t q);

void qf_mf4_convert(float4x4 m, quatf_t q);
void qd_md4_convert(double4x4 m, quatd_t q);
void qf_mf4_convert_inv(float4x4 m, quatf_t q);
void qd_md4_convert_inv(double4x4 m, quatd_t q);

quatf_t qf_slerp(quatf_t q0, quatf_t q1, float t);
quatd_t qd_slerp(quatd_t q0, quatd_t q1, double t);

/*!
 * \brief Converts a rotation matrix to a quaternion.
 * 
 * The m_q_convert function converts a rotational matrix m into a quaternion q.
 * The rotational matrix shall be specified for a right handed coordinate
 * system. For more information, see Ken Shoemake's paper on quaternion
 * rotation.
 * 
 * \param q A reference to the resulting quaternion.
 * \param m The rotational matrix.
 */

quatf_t mf4_qf_convert(float4x4 m) __attribute__ ((__nonnull__));
quatf_t mf3_qf_convert(float3x3 m) __attribute__ ((__nonnull__));
quatd_t md3_qd_convert(double3x3 m) __attribute__ ((__nonnull__));

quatf_t qf_add(const quatf_t a, const quatf_t b);
quatd_t qd_add(const quatd_t a, const quatd_t b);


quatf_t qf_mul(const quatf_t a, const quatf_t b);
quatd_t qd_mul(const quatd_t a, const quatd_t b);



quatf_t qf_s_mul(quatf_t q, float d);
quatd_t qd_s_mul(quatd_t q, double d);

quatf_t qf_s_div(const quatf_t q, float d);
quatd_t qd_s_div(const quatd_t q, double d);


float qf_dot(quatf_t a, quatf_t b) __attribute__ ((__pure__));
double qd_dot(quatd_t a, quatd_t b) __attribute__ ((__pure__));

    
float3 qf_cross(const quatf_t a, const quatf_t b);
double3 qd_cross(const quatd_t a, const quatd_t b);
    
float qf_abs(const quatf_t q);
double qd_abs(const quatd_t q);
    
quatf_t qf_conj(const quatf_t q);
quatd_t qd_conj(const quatd_t q);


quatf_t qf_repr(const quatf_t q);


quatf_t qf_div(const quatf_t a, const quatf_t b);
quatd_t qd_div(const quatd_t a, const quatd_t b);

/*!
 * \brief   Creates a rotation quaternion
 * 
 * The q_rot function creates a rotation quaternion for a right hand rotation
 * of alpha radians around the specified axis.
 * 
 * \param q     A reference to the new quaternion
 * \param axis  A unit vector describing the axis of rotation 
 * \param alpha Rotation in radians.
*/
quatf_t qf_rotv(float3 axis, float alpha);
quatd_t qd_rotv(double3 axis, double alpha);

quatf_t qf_rot(float x, float y, float z, float alpha);
quatd_t qd_rot(double x, double y, double z, double alpha);

#define Q_ROT_X(q, r)                                           \
    do {                                                        \
        float3 _v = {S_CONST(1.0), S_CONST(0.0),        \
                     S_CONST(0.0)};      \
        (q) = q_rotv(_v, r);                                    \
    } while (0)

#define Q_ROT_Y(q, r)                                           \
    do {                                                        \
        float3 _v = {S_CONST(0.0), S_CONST(1.0),        \
                    S_CONST(0.0)};      \
        (q) = q_rotv(_v, r);                                    \
    } while (0)


#define Q_ROT_Z(q, r)                                           \
    do {                                                        \
        float3 _v = {S_CONST(0.0), S_CONST(0.0),        \
                    S_CONST(1.0)};      \
        (q) = q_rotv(_v, r);                                    \
    } while (0)


quatf_t qf_normalise(quatf_t q);
quatd_t qd_normalise(quatd_t q);

static inline quatf_t
qf_vf3_rot(quatf_t q, float3 v, float dt)
{
  quatf_t qv = {v.x, v.y, v.z, 0.0f};
  quatf_t qp = qf_add(q, qf_s_mul(qf_mul(qv, q), dt/2.0));
  return qp;
}

static inline quatd_t
qd_vd3_rot(quatd_t q, double3 v, double dt)
{
  quatd_t qv = {v.x, v.y, v.z, 0.0};
  quatd_t qp = qd_add(q, qd_s_mul(qd_mul(qv, q), dt/2.0));
  return qp;
}


static inline quatf_t
qf_vd3_rot(quatf_t q, double3 v, double dt)
{
  quatf_t qv = {v.x, v.y, v.z, 0.0f};
  quatf_t qp = qf_add(q, qf_s_mul(qf_mul(qv, q), dt/2.0));
  return qp;
}


#ifdef __cplusplus
}
#endif 

#endif /* ! __QUATERNIONS_H__ */
