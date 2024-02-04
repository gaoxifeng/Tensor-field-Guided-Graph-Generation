#include "hierarchy.h"
#include "positions.h"
#include <math.h>
#include "timer.h"
void MultiResolutionHierarchy::smoothPositionsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic) {
    const SMatrix &L = mL[l];
    const MatrixXf &V = mV[l], &N = mN[l], &Q = mQ[l], &S = mS[l], &C = mC[l];
    MatrixXf &O = mO[l];

    Timer<> timer;

    double error = 0;
    int nLinks = 0;
    MatrixXf O_new(O.rows(), O.cols());
    tbb::spin_mutex mutex;

#if PARALLEL
    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) L.outerSize(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            std::vector<std::pair<uint32_t, Float>> neighbors;
            double errorLocal = 0;
            int nLinksLocal = 0;
            for (uint32_t k = range.begin(); k != range.end(); ++k) {
#endif				
			//std::vector<std::pair<uint32_t, Float>> neighbors;
			//double errorLocal = 0;
			//int nLinksLocal = 0;
			//for (uint32_t k = 0; k != L.outerSize(); ++k) {

                SMatrix::InnerIterator it(L, k);

                uint32_t i = it.row();

                const Vector3f q_i = Q.col(i);
                const Vector3f n_i = N.col(i);
                const Vector3f v_i = V.col(i);
                const Vector3f c_i = C.col(i);

                Vector3f o_i = O.col(i);

                if(l == 0 && BC[i])
                {
                    ;//O_new.col(i) = v_i;
                    //continue;
                }

                const Vector2f scale0 = S.col(i).head(2) * mScale;
                //Vector2f scale0(S(2, i) * mScale, S(2,i) * mScale);
                Vector2f Invscale0 (1.0/scale0[0], 1.0/scale0[1]);

                neighbors.clear();
                for (; it; ++it) {
                    uint32_t j = it.col();
                    if (i == j)
                        continue;
                    neighbors.push_back(std::make_pair(j, it.value()));
                }

                if (randomization && neighbors.size() > 0)
                    pcg32(mPositionIterations, k)
                        .shuffle(neighbors.begin(), neighbors.end());

                Float weightSum = 0.f;
                for (auto n : neighbors) {
                    uint32_t j = n.first;
                    Float value = n.second;

                    const Vector3f q_j = Q.col(j), v_j = V.col(j), n_j = N.col(j);
                    Vector3f o_j = O.col(j);
                    const Vector2f scale1 = S.col(j).head(2) * mScale;
                    //Vector2f scale1 (S(2, j) * mScale, S(2, j) * mScale);
                    Vector2f Invscale1 (1.0/scale1[0], 1.0/scale1[1]);

                    if (extrinsic) {
                        errorLocal += (O.col(i) -
                              PosyExtrinsic2D(O.col(i), q_i, n_i, v_i, scale0, Invscale0, o_j, q_j, n_j,
                                             v_j, scale1, Invscale1)).norm();
                        o_j = PosyExtrinsic2D(o_i, q_i, n_i, v_i, scale0, Invscale0, o_j, q_j, n_j,
                                             v_j, scale1, Invscale1);

                    } 
                    o_i = value * o_j + weightSum * o_i;
                    weightSum += value;
                    o_i /= weightSum;
                    nLinksLocal++;
                }

				
                o_i = PosyExtrinsic2D_round(o_i, q_i, n_i, v_i, scale0, Invscale0);
                o_i -= n_i.dot(o_i - v_i) * n_i;

				if (BV[l][i] && c_i != Vector3f::Zero())
                {
					//o_i = q_i.dot(o_i - v_i) * q_i + v_i;
                    o_i = c_i.dot(o_i - v_i) * c_i + v_i;
                    //cout<<"n_i: "<<n_i<<" c_i: "<<c_i <<" o_i: "<<o_i<<endl;
                }else
                    ;//cout<<"n_i: "<<n_i<<" o_i: "<<o_i<<endl;

				O_new.col(i) = o_i;
            }
#if PARALLEL
			tbb::spin_mutex::scoped_lock guard(mutex);
            error += errorLocal;
            nLinks += nLinksLocal;
			}
		);
#endif

    mOrientationIterations++;
	O = std::move(O_new);
}

void MultiResolutionHierarchy::smoothPositionsTet(uint32_t l, bool alignment, bool randomization) {
    const SMatrix &L = mL[l];
    const MatrixXf &V = mV[l], &N = mN[l], &Q = mQ[l], &C = mC[l], &S = mS[l];
    MatrixXf &O = mO[l];

    Timer<> timer;
    //timer.beginStage("Smoothing orientations at level " + std::to_string(l));

    double error = 0;
    int nLinks = 0;
    MatrixXf O_new(O.rows(), O.cols());
    tbb::spin_mutex mutex;

#if 1
    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) L.outerSize(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            std::vector<std::pair<uint32_t, Float>> neighbors;
            double errorLocal = 0;
            int nLinksLocal = 0;
            for (uint32_t k = range.begin(); k != range.end(); ++k) {
#endif
				// std::vector<std::pair<uint32_t, Float>> neighbors;
				// double errorLocal = 0;
				// int nLinksLocal = 0;
				// for (uint32_t k = 0; k <L.outerSize(); ++k) {
					
					SMatrix::InnerIterator it(L, k);

                uint32_t i = it.row();

                const Quaternion q_i = Q.col(i);
                const Vector3f n_i = N.col(i), v_i = V.col(i), c_i = C.col(i);
                Vector3f o_i = O.col(i);

                // if(l == 0 && BC[i])
                // {
                //     O_new.col(i) = v_i;
                //     continue;
                // }
                
                const Vector3f scale0 = S.col(i).head(3) * mScale;
                const Vector3f Invscale0 (1.0/scale0[0], 1.0/scale0[1], 1.0/scale0[2]);

                neighbors.clear();
                for (; it; ++it) {
                    uint32_t j = it.col();
                    if (i == j)
                        continue;
                    neighbors.push_back(std::make_pair(j, it.value()));
                }

                if (randomization && neighbors.size() > 0)
                    pcg32(mPositionIterations, k)
                        .shuffle(neighbors.begin(), neighbors.end());

                Float weightSum = 0.f;
                for (auto n : neighbors) {
                    uint32_t j = n.first;
                    Float value = n.second;

                    const Quaternion q_j = Q.col(j);
                    Vector3f o_j = O.col(j);

                    const Vector3f scale1 = S.col(j).head(3) * mScale;
                 
                    o_j = Posy3D(o_i, q_i,scale0, o_j, q_j, scale1, is_anisotropy).first;
                    o_i = value * o_j + weightSum * o_i;
                    weightSum += value;
                    o_i /= weightSum;

					if (alignment && n_i != Vector3f::Zero()) {
						double scale = scale0.maxCoeff(), invscale = 1.0 / scale;
						auto dp = n_i.dot(c_i - o_i) * invscale;
						o_i += (dp - round(dp)) * n_i * scale;
					}
				}
				o_i = Posy3DfindClosest(o_i, q_i, v_i, scale0, Invscale0);

                if(std::isnan(o_i[0]))
				{
                    cout<<n_i<< " "<<q_i<<endl;
                    cout<<o_i<< " "<<l<<" "<<i<<endl;

					cout << scale0 << " " << Invscale0 << " " << c_i << endl;

					for (auto n : neighbors) {
						uint32_t j = n.first;
						Float value = n.second;

						const Quaternion q_j = Q.col(j);
						Vector3f o_j = O.col(j);

						const Vector3f scale1 = S.col(j).head(3) * mScale;
						//const Vector3f scale0(S(3, i) * mScale, S(3,i) * mScale, S(3, i) * mScale);

						cout << " " << q_j << endl;
						cout << o_j << endl;

						cout << scale1  << endl;
					}


					std::cin.get(); 
                }
                
				O_new.col(i) = o_i;

			}
            tbb::spin_mutex::scoped_lock guard(mutex);

#if 1
        }
    );
#endif

    //timer.endStage("E = " + std::to_string(error / nLinks));
    mOrientationIterations++;
    O = std::move(O_new);
}


void MultiResolutionHierarchy::prolongPositions(int level) {
    
	const SMatrix &P = mP[level];
	for (int k = 0; k < P.outerSize(); ++k) {
		SMatrix::InnerIterator it(P, k);
		for (; it; ++it) {

			Vector3f o_i = mO[level + 1].col(it.col());;
			
			Vector3f v_i = mV[level].col(it.row());
			Vector3f n_i = mN[level].col(it.row());

			if (!D3) {
				Vector3f q_i = mQ[level].col(it.row()); 
				if (BV[level][it.row()]) {
					o_i = q_i.dot(o_i - v_i) * q_i + v_i;
				}
				o_i -= n_i.dot(o_i - v_i) * n_i;
				mO[level].col(it.row()) = o_i;
			}
			else {
				Vector3f s_i = mS[level].col(it.row());
				Vector3f c_i = mC[level].col(it.row());

				if (n_i != Vector3f::Zero()) {
					double scale = s_i.maxCoeff(), invscale = 1.0 / scale;
					auto dp = n_i.dot(c_i - o_i) * invscale;
					o_i += (dp - round(dp)) * n_i * scale;
				}

				mO[level].col(it.row()) = o_i;
			}
		}
	}
}
