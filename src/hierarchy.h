#pragma once

#include "common.h"
#include "aabb.h"
#include "lock.h"
#include "adjacency.h"
#include "meshio.h"
#include <nanogui/serializer/core.h>
#include <nanogui/serializer/sparse.h>
#include <unordered_map>
#include <algorithm> 
#include <queue>
#include <set>

using namespace std;

struct MeshStats {
	AABB mAABB;
	double mAverageEdgeLength = 0;
	double mMaximumEdgeLength  = 0;

	MeshStats() :
		mAverageEdgeLength(0.0f),
		mMaximumEdgeLength(0.0f){ }
    void compute_mesh_stats(const MatrixXu &E_, const MatrixXf &V_){
        for(int i=0;i<E_.cols();i++){
            Float edge_length = (V_.col(E_(0, i)) - V_.col(E_(1, i))).norm();
            mAverageEdgeLength += edge_length;
            mMaximumEdgeLength = std::max(mMaximumEdgeLength, (double)edge_length);
        }
        if(E_.cols()) mAverageEdgeLength/=E_.cols();
    }
};
class MultiResolutionHierarchy {
public:
    MultiResolutionHierarchy();

    bool load(const std::string &filename);
	//protected:
	void build();

    MatrixXf &V(uint32_t i = 0) { return mV[i]; }
    const MatrixXf &V(uint32_t i = 0) const { return mV[i]; }

    MatrixXf &N(uint32_t i = 0) { return mN[i]; }
    const MatrixXf &N(uint32_t i = 0) const { return mN[i]; }

    MatrixXf &Q(uint32_t i = 0) { return mQ[i]; }
    const MatrixXf &Q(uint32_t i = 0) const { return mQ[i]; }

    MatrixXf &O(uint32_t i = 0) { return mO[i]; }
    const MatrixXf &O(uint32_t i = 0) const { return mO[i]; }

    MatrixXf &C(uint32_t i = 0) { return mC[i]; }
    const MatrixXf &C(uint32_t i = 0) const { return mC[i]; }

    MatrixXu &E() { return mE; }
    const MatrixXu &E() const { return mE; }

    SMatrix &L(uint32_t i = 0) { return mL[i]; }
    const SMatrix &L(uint32_t i) const { return mL[i]; }

    void smoothOrientationsTet(uint32_t l, bool alignment, bool randomization);
    void smoothOrientationsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic);
    void prolongOrientations(int level);

    void smoothPositionsTet(uint32_t l, bool alignment, bool randomization);
    void smoothPositionsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic);
    
    void prolongPositions(int level);

    AABB aabb() const { return mAABB; }

    size_t vertexCount() const { return mV.size() > 0 ? mV[0].cols() : 0; }
    size_t faceCount() const { return mF.cols(); }
    size_t tetCount() const { return mT.cols(); }

    Float averageEdgeLength() const { return mAverageEdgeLength; }
	Float scale() const { return ratio_scale; }
	void setScale(Float scale) { 
		ratio_scale = scale; 
		mScale = diagonalLen * scale; 
		mInvScale = 1.f / mScale;
	}

    ordered_lock &mutex() const { return mMutex; }

    int levels() const { return mL.size(); }

    void construct_Graph(const Mesh &m, Graph &g);
	void construct_G2VE(const Graph &g_, MatrixXf &V);
	void construct_Es(const MatrixXu &mE, std::vector<tuple_E> &mEs);
    void composit_edges_colors(MatrixXf &Result_Vs, std::vector<tuple_E> &Es_to_render, MatrixXf &Result_edges);
	void composit_edges_colors(Graph &g_, MatrixXf &Result_edges);
    void edge_tagging2D();
    void edge_tagging3D();

    void graph_extraction();    
	void edge_tagging(Graph &g_);
	void detect_parallel_redundant_edges(Graph &g_, vector<int> &ledges);
	void split_parallel_redundant_edges(Graph &g_, vector<int> &ledges);
	void connect_close_points(Graph &g_);
    bool split_long_edges(Graph &g);
    bool collapse_fuse_edges(Graph &g);
    void remove_dangling(Graph &g_);

    void graph_extraction_fast();
    void group_representative(const Graph &g_, const vector<int> &vs, Graph_Node &gn, bool averaging = true);
    void split_group(Graph &g_, vector<vector<int>> & sets);
    void split_group_based_on_color(Graph &g_, vector<vector<int>> & sets);
    void grouping_collapse(Graph &g_);
	void re_coloring_diagonals(Graph &g_);
	void remove_diagonals(Graph &g_);
	void re_projection(Graph &g_);

    std::tuple<Edge_tag, Float, VectorXf> posy_info(Graph &g_, int v0, int v1, bool averaging = false);
	std::tuple<Edge_tag, Float, VectorXf> posy_info(int dim, Graph_Node &v0, Graph_Node &v1, int vn0, int vn1, bool averaging);
    std::tuple<Edge_tag, Float, VectorXf> posy_info(Graph &g_, int v0, int v1, Vector3f p0, Vector3f p1);
    void interpolate_twovs(VectorXf q0, Vector3f g0, Vector3f n0, Vector3f v0, Vector3f s0,
		VectorXf qj, Vector3f gj, Vector3f nj, Vector3f vj, Vector3f sj,
		VectorXf &qn, Vector3f &gn, Vector3f &nn, Vector3f &vn, Vector3f &sn,
		int dimension);

public:
    bool D3 = false;
	//for both 2D & 3D 
    std::vector<MatrixXf> mV;
    std::vector<MatrixXf> mN;
    std::vector<VectorXb> BV;
    std::vector<MatrixXf> mQ;
    std::vector<MatrixXf> mS;
    std::vector<MatrixXf> mO;
    std::vector<MatrixXf> mC;
    std::vector<SMatrix> mL;
    std::vector<SMatrix> mP;
    MatrixXu mE;
    MatrixXu mF;
    MatrixXu mT;
	bool Q_FROM_FILE = false;
	MatrixXf QoF;
    MatrixXf SoF;
    VectorXb BC;
    VectorXi BC_render;
    bool is_anisotropy = false;
    
    uint32_t mOrientationIterations;
    uint32_t mPositionIterations;
    AABB mAABB;
    Float mAverageEdgeLength;
    mutable ordered_lock mMutex;
    Float mScale, mInvScale;
	Float diagonalLen;
	Float ratio_scale;
public:
	std::string outpath;
	MeshStats ms;

    Mesh m;
	Mesh m_sur;
    Graph g;
    std::vector<tuple_E> mEs;

    int step=0;

	MatrixXf mV_final;
    MatrixXu mE_final;
	MatrixXf mQ_final;
	MatrixXf mN_final;
	//for rendering
	MatrixXf E_input_rend;
	MatrixXf mO_center, mV_tag_rend;
	MatrixXf mV_final_rend;
	MatrixXf E_rend, E_O_rend, E_I_rend;
	MatrixXf E_rend_o, E_O_rend_o, E_I_rend_o;
	MatrixXf E_tag_rend, E_tag_left_rend;
	MatrixXf E_final_rend;
};
