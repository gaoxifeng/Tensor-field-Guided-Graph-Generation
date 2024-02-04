#pragma once

#include "orientations.h"
#include <algorithm>

inline Vector3f middle_point(const Vector3f &p0, const Vector3f &n0, const Vector3f &p1, const Vector3f &n1) {
    Float n0p0 = n0.dot(p0), n0p1 = n0.dot(p1),
          n1p0 = n1.dot(p0), n1p1 = n1.dot(p1),
          n0n1 = n0.dot(n1),
          denom = 1.0f / (1.0f - n0n1*n0n1 + 1e-4f),
          lambda_0 = 2.0f*(n0p1 - n0p0 - n0n1*(n1p0 - n1p1))*denom,
          lambda_1 = 2.0f*(n1p0 - n1p1 - n0n1*(n0p1 - n0p0))*denom;
    return 0.5f * (p0 + p1) - 0.25f * (n0 * lambda_0 + n1 * lambda_1);
}

inline Vector3f
PosyExtrinsic2D(const Vector3f &o0_, const Vector3f &q0, const Vector3f &n0, const Vector3f &p0,
						 const Vector2f scale0, const Vector2f Invscale0,
                         const Vector3f &o1_, const Vector3f &q1, const Vector3f &n1, const Vector3f &p1,
						 const Vector2f scale1, const Vector2f Invscale1) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;
    Matrix M0 = (Matrix() << q0, n0.cross(q0)).finished();
    Matrix M1 = (Matrix() << q1, n1.cross(q1)).finished();
    Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector2f dim0 = M0.transpose() * (middle - o0_);
	dim0[0] =dim0[0]*Invscale0[0];
	dim0[1] =dim0[1]*Invscale0[1];
	dim0 = floor (dim0);
	dim0[0] =dim0[0]*scale0[0];
	dim0[1] =dim0[1]*scale0[1];
	Vector3f rel0_ =  M0 * dim0;
	
	Vector2f dim1 = M1.transpose() * (middle - o1_);
	dim1[0] =dim1[0]*Invscale1[0];
	dim1[1] =dim1[1]*Invscale1[1];
	dim1 = floor (dim1);
	dim1[0] =dim1[0]*scale1[0];
	dim1[1] =dim1[1]*scale1[1];
	Vector3f rel1_ =  M1 * dim1;

    Float best_score = std::numeric_limits<Float>::infinity();
    Vector3f o1 = o1_;

    for (int i=0; i<16; ++i) {
		Vector2f scale(((i & (1 << 0)) >> 0) * scale0[0], ((i & (1 << 1)) >> 1) * scale0[1]);

        Vector3f rel0 = rel0_ + M0 * scale;

		scale[0] = ((i & (1 << 2)) >> 2) * scale1[0];
		scale[1] = ((i & (1 << 3)) >> 3) * scale1[1];

        Vector3f rel1 = rel1_ + M1 * scale;

        Float score = (o0_ + rel0 - (o1_ + rel1)).squaredNorm();

        if (score < best_score) {
            best_score = score;
            o1 = o1_ + rel1 - rel0;
        }
    }

    return o1;
}
inline Vector3f PosyExtrinsic2D_round(const Vector3f &o, const Vector3f &q, const Vector3f &n,
                            const Vector3f &ref, const Vector2f scale0, const Vector2f Invscale0) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;
    Matrix M = (Matrix() << q, n.cross(q)).finished();
    Vector2f rel = M.transpose() * (ref - o);
	rel[0] = rel[0] * Invscale0[0];
	rel[1] = rel[1] * Invscale0[1];
	rel = round(rel);
	rel[0] = rel[0] * scale0[0];
	rel[1] = rel[1] * scale0[1];
	
    return o + M * rel;
}
inline std::tuple<Edge_tag, Float, Vector2f> posy2D_completeInfo(const Vector3f &o0_, const Vector3f &q0, const Vector3f &n0, const Vector3f &p0,
						 const Vector2f scale0,
                         const Vector3f &o1_, const Vector3f &q1, const Vector3f &n1, const Vector3f &p1,
						 const Vector2f scale1, bool anisotropy = false) {
	typedef Eigen::Matrix<Float, 3, 2> Matrix;


	Vector3f qn = (q0 + applyRotation(q0, n0, q1, n1)).normalized();
	
	auto scale1_ = scale1;
	if(anisotropy)
	{
		int which = findRotation(q0, n0, q1, n1);
		if (which == 1 || which == 3)
			std::swap(scale1_[0], scale1_[1]);
	}

	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o1 = rotateVectorIntoPlane(o1_ - middle, n1, n0) + middle;

	Matrix M = (Matrix() << qn, n0.cross(qn)).finished();

	Vector2f aveScale = (scale0 + scale1_)* 0.5;
	Vector2f aveInvScale(1.0/aveScale[0],1.0/aveScale[1]);

	Vector2f vv = M.transpose() * (o0_ - o1);
	vv[0] = vv[0] * aveInvScale[0];
	vv[1] = vv[1] * aveInvScale[1];

	for (uint32_t i = 0; i < 2; i++) vv[i] = std::abs(vv[i]);

	Edge_tag e_color;
	int res = 0; Float Weight=0; std::vector<bool> bits(2);
	for (uint32_t i = 0; i < 2; i++) {
		if (vv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 2; i++) Weight += vv[i] * vv[i];
		e_color = Edge_tag::R;
	}
	else if (res == 1) {
		for (int i = 0; i < 2; i++) {
			if (bits[i]) { Weight += (1 - vv[i]) * (1 - vv[i]); }
			else Weight += vv[i] * vv[i];
		}
		Weight += 1.0;
		e_color = Edge_tag::B;
	}
	else if (res == 2) {
		for (int i = 0; i < 2; i++)
			Weight += (1 - vv[i]) * (1 - vv[i]);
		Weight += std::sqrt(2.0);
		e_color = Edge_tag::D;
	}

	return std::make_tuple(e_color, Weight, vv);
}

inline std::pair<Vector3f, Vector3i>
Posy3D(const Vector3f &o0, const Quaternion &q0,
	   const Vector3f scale0,
       const Vector3f &o1, const Quaternion &q1,
       const Vector3f scale1, bool anisotropy = false) {

	Quaternion q1_ =Quaternion::applyRotation(q1, q0);
    Quaternion qn = (q0 + q1_).normalized();
    Matrix3f M = qn.toMatrix();

	Vector3f scale = scale1;

	if(anisotropy)
	{
		Eigen::Matrix<Float, 3, 3> M0, M1;
		M0 = q1.toMatrix();
		M1 = q1_.toMatrix();
		for(int j=0;j<3;j++){
			int id=0;Float dmax=-1;
			for(int k=0;k<3;k++){
				Float d = std::abs(M0.col(j).dot(M1.col(k)));
				if(d>dmax){dmax = d; id = k;}	
			}
			scale[j] = scale1[id];
		}
	}

	Vector3f aveScale = (scale0 + scale)* 0.5;
	Vector3f aveInvScale(1.0/aveScale[0],1.0/aveScale[1],1.0/aveScale[2]);
	Vector3f vv = M.transpose() * (o0 - o1);
	vv[0] = vv[0] * aveInvScale[0];
	vv[1] = vv[1] * aveInvScale[1];
	vv[2] = vv[2] * aveInvScale[2];
	Vector3f rel = round(vv);
	
	//vv = M * rel;
	vv = rel;
	vv[0] = vv[0] * aveScale[0];
	vv[1] = vv[1] * aveScale[1];
	vv[2] = vv[2] * aveScale[2];
	
	//return std::make_pair(o1 + vv, rel.cast<int>());
	return std::make_pair(o1 + M * vv, rel.cast<int>());
}

inline Vector3f Posy3DfindClosest(const Vector3f &o, const Quaternion &q, 
                            const Vector3f &ref, const Vector3f scale, const Vector3f Invscale) {
    Matrix3f M = q.toMatrix();
	
	Vector3f vv = M.transpose() * (ref - o);
	vv[0] = vv[0] * Invscale[0];
	vv[1] = vv[1] * Invscale[1];
	vv[2] = vv[2] * Invscale[2];
	Vector3f rel = round(vv);

	//vv = M * rel;
	vv = rel;
	vv[0] = vv[0] * scale[0];
	vv[1] = vv[1] * scale[1];
	vv[2] = vv[2] * scale[2];

    return o + M * vv;
}
inline std::tuple<Edge_tag, Float, Vector3f> posy3D_completeInfo(const Vector3f &o0, const Quaternion &q0, const Vector3f scale0, 
	const Vector3f &o1, const Quaternion &q1, const Vector3f scale1, bool anisotropy = false) {

	Quaternion q1_ =Quaternion::applyRotation(q1, q0);
	Quaternion qn = (q0 + q1_).normalized();
	Matrix3f M = qn.toMatrix();

	Vector3f scale = scale1;;
	if(anisotropy)
	{
		Eigen::Matrix<Float, 3, 3> M0, M1;
		M0 = q1.toMatrix();
		M1 = q1_.toMatrix();
		Vector3f scale;
		for(int j=0;j<3;j++){
			int id=0;Float dmax=-1;
			for(int k=0;k<3;k++){
				Float d = std::abs(M0.col(j).dot(M1.col(k)));
				if(d>dmax){dmax = d; id = k;}	
			}
			scale[j] = scale1[id];
		}
	}

	Vector3f aveScale = (scale0 + scale)* 0.5;
	Vector3f aveInvScale(1.0/aveScale[0],1.0/aveScale[1],1.0/aveScale[2]);
	
	Vector3f vvv = M.transpose() * (o0 - o1);
	vvv[0] = vvv[0] * aveInvScale[0];
	vvv[1] = vvv[1] * aveInvScale[1];
	vvv[2] = vvv[2] * aveInvScale[2];
	for (uint32_t i = 0; i < 3; i++) {
		vvv[i] = std::abs(vvv[i]);
	}

	Edge_tag e_color;
	short res = 0;  Float Weight = 0; std::vector<bool> bits(3);
	for (uint32_t i = 0; i < 3; i++) {
		if (vvv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 3; i++) Weight += vvv[i] * vvv[i];
		e_color = Edge_tag::R;
	}
	else if (res == 1) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 1.0;
		e_color = Edge_tag::B;
	}
	else if (res == 2) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 2.0;
		e_color = Edge_tag::D;
	}
	else if (res == 3) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 3.0;
		e_color = Edge_tag::H;
	}

	return std::make_tuple(e_color, Weight, vvv);
}
//==============================================================================
inline std::pair<Vector3f, Vector3i>
findClosestPair(const Vector3f &o0, const Quaternion &q0,
                const Vector3f &o1, const Quaternion &q1,
                Float scale, Float invScale) {
    Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
    Matrix3f M = qn.toMatrix();
    Vector3f rel = round(Vector3f(M.transpose() * ((o0 - o1) * invScale)));
    return std::make_pair(o1 + (M * rel) * scale, rel.cast<int>());
}
inline  Vector3f exact_3dir(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	return M.transpose() * ((o0 - o1) * invScale);
	
}
inline std::tuple<short, Float, Vector3f> posy3D_completeInfo(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	Vector3f vvv = (M.transpose() * ((o0 - o1) * invScale));

	// we want to count how many coords of vvv will be rounded to 1
	for (uint32_t i = 0; i < 3; i++) {
		vvv[i] = std::abs(vvv[i]);
	}
	short res = 0;  Float Weight = 0; std::vector<bool> bits(3);
	for (uint32_t i = 0; i < 3; i++) {
		if (vvv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 3; i++) Weight += vvv[i] * vvv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 2.0;
	}
	else if (res == 3) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 3.0;
	}

	return std::make_tuple(res, Weight, vvv);
}

inline std::pair<short, short> assignColorBiased(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	Vector3f vvv = (M.transpose() * ((o0 - o1) * invScale));

	// we want to count how many coords of vvv will be rounded to 1
	for (uint32_t i = 0; i < 3; i++) {
		vvv[i] = std::abs(vvv[i]);
	}
	short res = 0, ind; std::vector<bool> bits(3);
	for (uint32_t i = 0; i < 3; i++) {
		if (vvv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		ind=-1;// for (int i = 0; i < 3; i++) if (vvv[i] > 0.2) res = 1;
	}
	else if (res == 1) {
		for (short i = 0; i < 3; i++) if (bits[i]) ind = i;
	}
	else if (res == 2) {
		for (short i = 0; i < 3; i++) if (!bits[i]) ind = i;
	}else if(res==3)
		ind = -1;
	return std::make_pair(res,ind);
}
inline std::pair<short, Float> assignColorWeighted3D(const Vector3f &o0, const Quaternion &q0,
	const Vector3f &o1, const Quaternion &q1,
	Float scale, Float invScale) {
	Quaternion qn = (q0 + Quaternion::applyRotation(q1, q0)).normalized();
	Matrix3f M = qn.toMatrix();
	Vector3f vvv = (M.transpose() * ((o0 - o1) * invScale));

	// we want to count how many coords of vvv will be rounded to 1
	for (uint32_t i = 0; i < 3; i++) {
		vvv[i] = std::abs(vvv[i]);
	}
	short res = 0;  Float Weight = 0; std::vector<bool> bits(3);
	for (uint32_t i = 0; i < 3; i++) {
		if (vvv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 3; i++) Weight += vvv[i] * vvv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 2.0;
	}
	else if (res == 3) {
		for (int i = 0; i < 3; i++) {
			if (bits[i]) { Weight += (1 - vvv[i]) * (1 - vvv[i]); }
			else Weight += vvv[i] * vvv[i];
		}
		Weight += 3.0;
	}

	return std::make_pair(res, Weight);
}


inline Vector3f
findClosestPairExtrinsic(const Vector3f &o0_, const Vector3f &q0, const Vector3f &n0, const Vector3f &p0,
                         const Vector3f &o1_, const Vector3f &q1, const Vector3f &n1, const Vector3f &p1,
                         Float scale, Float invScale) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;
    Matrix M0 = (Matrix() << q0, n0.cross(q0)).finished();
    Matrix M1 = (Matrix() << q1, n1.cross(q1)).finished();
    Vector3f middle = middle_point(p0, n0, p1, n1);
    Vector3f rel0_ =  M0 * floor(Vector2f(M0.transpose() * (middle - o0_) * invScale)) * scale;
    Vector3f rel1_ =  M1 * floor(Vector2f(M1.transpose() * (middle - o1_) * invScale)) * scale;

    Float best_score = std::numeric_limits<Float>::infinity();
    Vector3f o1 = o1_;

    for (int i=0; i<16; ++i) {
        Vector3f rel0 = rel0_ + M0 * Vector2f((i & (1 << 0)) >> 0,
                                              (i & (1 << 1)) >> 1) * scale;
        Vector3f rel1 = rel1_ + M1 * Vector2f((i & (1 << 2)) >> 2,
                                              (i & (1 << 3)) >> 3) * scale;

        Float score = (o0_ + rel0 - (o1_ + rel1)).squaredNorm();

        if (score < best_score) {
            best_score = score;
            o1 = o1_ + rel1 - rel0;
        }
    }

    return o1;
}

inline std::pair<Vector3f, Vector2i>
findClosestPair(const Vector3f &o0,  const Vector3f &q0,  const Vector3f &n0, const Vector3f &p0,
                const Vector3f &o1_, const Vector3f &q1_, const Vector3f &n1, const Vector3f &p1,
                Float scale, Float invScale) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;

    Vector3f qn = (q0 + applyRotation(q0, n0, q1_, n1)).normalized();
    Vector3f middle = middle_point(p0, n0, p1, n1);
    Vector3f o1 = rotateVectorIntoPlane(o1_ - middle, n1, n0) + middle;

    Matrix M = (Matrix() << qn, n0.cross(qn)).finished();
    Vector2f rel = round(Vector2f(M.transpose() * ((o0 - o1) * invScale)));
    return std::make_pair(o1 + (M * rel) * scale, rel.cast<int>());
}
inline std::pair<int, Float> assignColorWeighted2D(const Vector3f &o0, const Vector3f &q0, const Vector3f &n0, const Vector3f &p0,
	const Vector3f &o1_, const Vector3f &q1_, const Vector3f &n1, const Vector3f &p1,
	Float scale, Float invScale) {
	typedef Eigen::Matrix<Float, 3, 2> Matrix;

	Vector3f qn = (q0 + applyRotation(q0, n0, q1_, n1)).normalized();
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o1 = rotateVectorIntoPlane(o1_ - middle, n1, n0) + middle;

	Matrix M = (Matrix() << qn, n0.cross(qn)).finished();
	//Vector2f rel = round();
	Vector2f vv = Vector2f(M.transpose() * ((o0 - o1) * invScale));
	for (uint32_t i = 0; i < 2; i++) vv[i] = std::abs(vv[i]);

	short res = 0; Float Weight=0; std::vector<bool> bits(2);
	for (uint32_t i = 0; i < 2; i++) {
		if (vv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 2; i++) Weight += vv[i] * vv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 2; i++) {
			if (bits[i]) { Weight += (1 - vv[i]) * (1 - vv[i]); }
			else Weight += vv[i] * vv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 2; i++)
			Weight += (1 - vv[i]) * (1 - vv[i]);
		Weight += std::sqrt(2.0);
	}

	return std::make_pair(res, Weight);
}

inline std::tuple<int, Float, Vector2f> posy2D_completeInfo_global(const Vector3f &o0, const Vector3f &o1, Float scale, Float invScale) {
	typedef Eigen::Matrix<Float, 3, 2> Matrix;
	Vector2f vv = ((round(o0) - round(o1)) * 1).segment(0, 2);
	for (uint32_t i = 0; i < 2; i++) vv[i] = std::abs(vv[i]);

	short res = 0; Float Weight = 0; std::vector<bool> bits(2);
	for (uint32_t i = 0; i < 2; i++) {
		if (vv[i] > 0.5) { res++; bits[i] = true; }
		else bits[i] = false;
	}
	if (res == 0) { // it's a red! 
		for (int i = 0; i < 2; i++) Weight += vv[i] * vv[i];
	}
	else if (res == 1) {
		for (int i = 0; i < 2; i++) {
			if (bits[i]) { Weight += (1 - vv[i]) * (1 - vv[i]); }
			else Weight += vv[i] * vv[i];
		}
		Weight += 1.0;
	}
	else if (res == 2) {
		for (int i = 0; i < 2; i++)
			Weight += (1 - vv[i]) * (1 - vv[i]);
		Weight += std::sqrt(2.0);
	}

	return std::make_tuple(res, Weight, vv);
}
inline Vector3f findClosest(const Vector3f &o, const Quaternion &q, 
                            const Vector3f &ref, Float scale, Float invScale) {
    Matrix3f M = q.toMatrix();
    Vector3f rel = M.transpose() * ((ref - o) * invScale);
    return o + (M * round(rel)) * scale;
}

inline Vector3f findClosest(const Vector3f &o, const Vector3f &q, const Vector3f &n,
                            const Vector3f &ref, Float scale, Float invScale) {
    typedef Eigen::Matrix<Float, 3, 2> Matrix;
    Matrix M = (Matrix() << q, n.cross(q)).finished();
    Vector2f rel = M.transpose() * ((ref - o) * invScale);
    return o + (M * round(rel)) * scale;
}


//projection
void projectPointOnQuad(const vector<Vector3d>& quad_vs, vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN);
void projectPointOnTriangle(const vector<Vector3d>& tri_vs, const vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN);
template <typename T>
T bilinear(const T& v1, const T& v2, const T& v3, const T& v4, const Vector2d& uv);
template <typename T>
T barycentric(const T& v1, const T& v2, const T& v3, const Vector2d& uv);
template <typename T>
bool num_equal(const T& x, const T& y, const double &precision);