#include "hierarchy.h"
#include "meshio.h"
#include "timer.h"
#include "quat.h"
#include "bvh.h"
#include "orientations.h"

MultiResolutionHierarchy::MultiResolutionHierarchy() {
    mV = { MatrixXf::Zero(3, 1) };
    mN = { MatrixXf::Zero(3, 1) };
    mO = { MatrixXf::Zero(3, 1) };
    mQ = { MatrixXf::Zero(4, 1) };
	ratio_scale = 3.0;
}

bool MultiResolutionHierarchy::load(const std::string &filename) {

    mV.resize(1);
	mV[0] = MatrixXf::Zero(3, 1);
	BV.resize(1);
	mN.resize(1);
	mC.resize(1);
	
	if(mesh_preprocessing(filename, mE, mV[0], BV[0], mN[0], mC[0], QoF, SoF, BC, m)){	
		Q_FROM_FILE = true;
	}else return false;

	//write_surface_mesh_OBJ(m, filename + "_.obj");
	D3 = m.type == Mesh_type::Tri? false:true;

	Float min = SoF.block(0,0, 2, SoF.cols()).minCoeff(), max = SoF.block(0,0, 2, SoF.cols()).maxCoeff();	
	//for (int i = 0; i < SoF.rows(); i++)
	//	for (int j = 0; j < SoF.cols(); j++)
	//		SoF(i, j) = 1.0 / SoF(i, j) * max;

	if (D3)
	{
		extract_surface_mesh(m, m_sur);
		vector<vector<uint32_t>> HF_(m_sur.Fs.size());
		for (auto &f : m_sur.Fs)
			HF_[f.id] = f.vs;
		orient_polygon_mesh(m_sur.V, HF_);
		//write_surface_mesh_OBJ(m_sur, path + "sur.obj");
	}


	cout<<"max, min: "<<max<<" "<<min<<endl;
	is_anisotropy = false;
	if(min == max)
		is_anisotropy = true;

	mV.resize(1);
	mAABB = AABB(
		mV[0].rowwise().minCoeff(),
		mV[0].rowwise().maxCoeff()
	);

	ms.compute_mesh_stats(mE, mV[0]);
	diagonalLen = (mAABB.max - mAABB.min).norm();
	ratio_scale = ms.mAverageEdgeLength * 3.0 / diagonalLen;

	construct_Es(mE, mEs);
	composit_edges_colors(mV[0], mEs, E_input_rend);

	BC_render.resize(BC.size());
	for(int i=0;i<BC.size();i++)
		BC_render[i] = BC[i];
	std::cout<<"boundary contraints: "<<BC_render.size()<<" "<<BC_render.sum()<<std::endl;
    return true;
}
void MultiResolutionHierarchy::build() {

	Timer<> timer;
	mV.resize(1);
	BV.resize(1);	
	mN.resize(1);
	mC.resize(1);

	struct WeightedEdge {
		WeightedEdge(uint32_t _i0, uint32_t _i1, Float weight)
			: weight(weight), i0(_i0), i1(_i1) {
			if (i0 > i1)
				std::swap(i0, i1);
		}

		bool operator<(const WeightedEdge &e) const {
			return std::tie(weight, i0, i1) < std::tie(e.weight, e.i0, e.i1);
		}

		Float weight;
		uint32_t i0, i1;
	};
	//laplacian matrix
	mL.clear(); mP.clear();

	std::vector<Triplet> triplets;
	for (int i = 0; i<mE.cols();i++){
		triplets.push_back(Triplet(mE(0, i), mE(1, i), 1.f));
		triplets.push_back(Triplet(mE(1, i), mE(0, i), 1.f));
	}
	mL.resize(1);
	mL[0].resize(mV[0].cols(), mV[0].cols());
	mL[0].setFromTriplets(triplets.begin(), triplets.end());

	for (uint32_t i = 0; i < (uint32_t)mL[0].rows(); ++i) {
		Float sum = 1 / mL[0].row(i).sum();
		mL[0].row(i) *= sum;
		mL[0].coeffRef(i, i) = -sum;
	}
	mL[0].makeCompressed();

	timer.beginStage("Building hierarchy");
	std::vector<SMatrix> mR;
	while (mL[mL.size() - 1].cols() > 1) {
		const MatrixXf &V = mV[mV.size() - 1];
		const MatrixXf &N = mN[mN.size() - 1];
		const MatrixXf &C = mC[mC.size() - 1];
		const VectorXb &VB = BV[BV.size() - 1];
		const SMatrix &L = mL[mL.size() - 1];
		std::vector<bool> collapsed(L.cols(), false);
		std::vector<bool> visited(L.cols(), false);
		std::set<WeightedEdge> edges;

		double edgeSum = 0;
		size_t edgeCount = 0;
		for (int k = 0; k < L.outerSize(); ++k) {
			for (SMatrix::InnerIterator it(L, k); it; ++it) {
				if (it.col() == it.row())
					continue;
				Float length = (V.col(it.row()) - V.col(it.col())).norm();
				edgeSum += length;
				edgeCount += 1;
				edges.insert(WeightedEdge(it.row(), it.col(), length));
			}
		}
		if (mL.size() == 1)
			mAverageEdgeLength = edgeSum / edgeCount;

		std::vector<Triplet> P_triplets, R_triplets;
		std::vector<Vector3f> V_next, N_next, C_next;
		std::map<uint32_t, uint32_t> vertex_map;

		uint32_t nVertices = 0; VectorXb vb_flag(V.cols());vb_flag.setZero();
		for (auto const &e : edges) {
			visited[e.i0] = visited[e.i1] = true;
			if (collapsed[e.i0] || collapsed[e.i1])
				continue;
			collapsed[e.i0] = true;
			collapsed[e.i1] = true;
			P_triplets.push_back(Triplet(e.i0, nVertices, 1.0f));
			P_triplets.push_back(Triplet(e.i1, nVertices, 1.0f));
			R_triplets.push_back(Triplet(nVertices, e.i0, 0.5f));
			R_triplets.push_back(Triplet(nVertices, e.i1, 0.5f));
			V_next.push_back(0.5f * (V.col(e.i0) + V.col(e.i1)));

			if (VB[e.i0] || VB[e.i1]) vb_flag[nVertices] = true;

			Vector3f n = N.col(e.i0) + N.col(e.i1);
			Vector3f c = C.col(e.i0) + C.col(e.i1);
			if (N.col(e.i0) != Vector3f::Zero() &&
				N.col(e.i1) != Vector3f::Zero()) {
				if (D3) {
					n = N.col(e.i0).normalized();
					c = C.col(e.i0);
				}
				else {
					n.normalize();
					//c *= 0.5f;
					if(c != Vector3f::Zero())
						c.normalize();
				}
			}

			if(std::isnan(n[0]))
			{
				cout<<"error Here"<<mL.size()<<endl;
			}

			N_next.push_back(n);
			C_next.push_back(c);

			vertex_map[e.i0] = nVertices;
			vertex_map[e.i1] = nVertices;
			nVertices++;
		}

		for (uint32_t i = 0; i<V.cols(); ++i) {
			if (collapsed[i] || !visited[i])
				continue;
			P_triplets.push_back(Triplet(i, nVertices, 1.0f));
			R_triplets.push_back(Triplet(nVertices, i, 1.0f));
			V_next.push_back(V.col(i));
			N_next.push_back(N.col(i));
			C_next.push_back(C.col(i));
			vertex_map[i] = nVertices;

			if (VB[i]) vb_flag[nVertices] = true;

			nVertices++;
		}
		vb_flag.resize(nVertices);

		if (mL.size() != 1)
			std::cout << ", ";
		std::cout << nVertices;
		std::cout.flush();

		SMatrix P(V.cols(), nVertices), R(nVertices, V.cols());

		P.setFromTriplets(P_triplets.begin(), P_triplets.end());
		R.setFromTriplets(R_triplets.begin(), R_triplets.end());

		SMatrix L2 = R*L*P;
		MatrixXf V2(3, nVertices), N2(3, nVertices), C2(3, nVertices), Q2(4, nVertices);
		for (uint32_t i = 0; i<nVertices; ++i) {
			V2.col(i) = V_next[i];
			N2.col(i) = N_next[i];
			C2.col(i) = C_next[i];
		}

		BV.push_back(vb_flag);
		mR.push_back(std::move(R));
		mP.push_back(std::move(P));
		mN.push_back(std::move(N2));
		mV.push_back(std::move(V2));
		mC.push_back(std::move(C2));
		mL.push_back(L2);
	}
	std::cout << " ";
	timer.endStage();

	mQ.resize(mL.size());
	mS.resize(mL.size());
	mO.resize(mL.size());

	if (Q_FROM_FILE){
		mQ[0] = QoF;
		mS[0] = SoF;
	}

	pcg32 rng;

	for (uint32_t i = 0; i < mL.size(); ++i) {
		mO[i].resize(3, mV[i].cols());

		if(D3){
			for (uint32_t j = 0; j < mV[i].cols(); ++j)
			{
				//mO[i].col(j) = aabbRand(mAABB, rng);
				mO[i].col(j) = mV[i].col(j);
			}
		}else{
			for (uint32_t j = 0; j < mV[i].cols(); ++j) {
				// Vector3f n = mN[i].col(j), v = mV[i].col(j);
				// rng.nextFloat();
				// Vector3f o = aabbRand(mAABB, rng);
				// o -= n.dot(o - v) * n;
				// if(i == 0 && BC[j])
				// 	mO[i].col(j) = mV[i].col(j);
				// else
				// 	mO[i].col(j) = o;
				mO[i].col(j) = mV[i].col(j);
			}
		}

		if (Q_FROM_FILE && i == 0)continue;
		else {
			if(D3)
				mQ[i].resize(4, mV[i].cols());
			else 
				mQ[i].resize(3, mV[i].cols());
			mS[i].resize(3, mV[i].cols());
		}
	}
	//propagate up
	cout << "mQ " << mQ[0].cols() << "; mV" << mV[0].cols() << endl;
	for (uint32_t i = 1; i < mL.size(); ++i) {
		if (Q_FROM_FILE) {
			for (int k = 0; k < mR[i - 1].outerSize(); ++k) {
				SMatrix::InnerIterator it(mR[i - 1], k);
				mQ[i].col(it.row()).setZero();
				mS[i].col(it.row()).setZero();
				Vector3f n_i = mN[i].col(it.row());

				int in = 0;				

				if(D3){
					int fist = -1;
					for (; it; ++it) {
						Quaternion q_j = mQ[i - 1].col(it.col());

						Vector3f scale = mS[i - 1].col(it.col());
						if(in!=0){

    	    	            Quaternion q_j_ = Quaternion::applyRotation(q_j, mQ[i-1].col(fist));
							Eigen::Matrix<Float, 3, 3> M0, M1;
							M0 = q_j.toMatrix();
							M1 = q_j_.toMatrix();
							for(int j=0;j<3;j++){
								int id=0;Float dmax=-1;
								for(int m=0;m<3;m++){
									Float d = std::abs(M0.col(j).dot(M1.col(m)));
									if(d>dmax){dmax = d; id = m;}	
								}
								scale[j] = mS[i - 1].col(it.col())[id];
							}

							mQ[i].col(it.row()) += q_j_; 
						}else {
							mQ[i].col(it.row()) += q_j;
							fist = it.col();
						}
						mS[i].col(it.row()) += scale;
						in++;
					}
					if(in != 0) mS[i].col(it.row()) /=in;
					else cout<<"error "<<endl;
					
					if (mQ[i].col(it.row()) == Eigen::Vector4f::Zero()) {
						mQ[i].col(it.row()) = Quaternion::Random(rng);
					}else mQ[i].col(it.row()).normalize();
					
				}else {	
					for (; it; ++it) {
						Vector3f q_j = mQ[i - 1].col(it.col());
						Vector3f n_j = mN[i - 1].col(it.col());

						int which = 0;
						if(in!=0){
							which = findRotation(mQ[i].col(it.row()), n_i, q_j, n_j);
							Vector3f q_j_ = applyRotationKeep(mQ[i].col(it.row()), n_i, q_j, n_j);
							mQ[i].col(it.row()) += q_j_ - n_i * n_i.dot(q_j_); 
						}else  mQ[i].col(it.row()) += q_j - n_i * n_i.dot(q_j);

						if(which==1||which ==3){
							mS[i].col(it.row())[0] += mS[i - 1].col(it.col())[1];
							mS[i].col(it.row())[1] += mS[i - 1].col(it.col())[0];
							mS[i].col(it.row())[2] += mS[i - 1].col(it.col())[2];
						}else mS[i].col(it.row()) += mS[i - 1].col(it.col());
						in++;
					}
					if(in != 0) mS[i].col(it.row()) /=in;
					else cout<<"error "<<endl;
					//if(mS[i].col(it.row())[2] == 0) cout<<"confusing"<<endl;

					if (mQ[i].col(it.row()) == Eigen::Vector3f::Zero()) {
						Vector3f n = mN[i].col(it.row()), v = mV[i].col(it.row());
						Vector3f s, t;
						coordinate_system(n, s, t);
						float angle = rng.nextFloat() * 2 * M_PI;
						mQ[i].col(it.row()) = s * std::cos(angle) + t * std::sin(angle);
					}else mQ[i].col(it.row()).normalize();

				}

			}
		}
	}

	mOrientationIterations = 0;
	mPositionIterations = 0;

 	mScale = diagonalLen * ratio_scale;
	//mScale = 5;
	mInvScale = 1.f / mScale;
	cout<<"mScale" << mScale <<endl;
	cout<<"mInvScale" << mInvScale <<endl;
}