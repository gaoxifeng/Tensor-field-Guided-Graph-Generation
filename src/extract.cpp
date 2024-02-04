#include "hierarchy.h"
#include "positions.h"
#include "bvh.h"
class BVH;
//2D&3D========================================================================================================//
std::priority_queue<tuple_E, std::vector<tuple_E>, LessThan> Es_reddash;
std::vector<uint32_t> V_map;
std::vector<std::vector<uint32_t>> Reverse_V_map;

std::vector<int> mV_flag;
std::vector<bool> mE_flag;

//2D===========================================================================================================//
void MultiResolutionHierarchy::construct_Graph(const Mesh &m_, Graph &g_)
{
	g_.dimension = m_.type == Mesh_type::Tri? 2:3;
	//node
	for(auto &v: m_.Vs)
	{
		Graph_Node gn;
		gn.id = g_.Vs.size();
		gn.boundary = v.boundary;
		gn.scale = v.scale;
		gn.v = m_.V.col(v.id);
		gn.position = mO[0].col(v.id);
		gn.direction = v.direction;
		gn.normal = v.normal;
		g_.Vs.push_back(gn);
	}
	//edge
	for(auto &e:m_.Es)
	{
		Graph_Edge ge;
		ge.vs.push_back(e.vs[0]);
		ge.vs.push_back(e.vs[1]);
		ge.boundary = e.boundary;
		ge.id = g_.Es.size();
		ge.time = 0;
		g_.Es.push_back(ge);

		g_.Vs[ge.vs[0]].nvs.insert(ge.vs[1]);
		g_.Vs[ge.vs[1]].nvs.insert(ge.vs[0]);
		g_.Vs[ge.vs[0]].nes.insert(ge.id);
		g_.Vs[ge.vs[1]].nes.insert(ge.id);
	}
}
void MultiResolutionHierarchy::construct_G2VE(const Graph &g_, MatrixXf &V) 
{
	V.resize(3, g_.Vs.size());
	for (auto &v : g_.Vs)
	{
		V.col(v.id) = v.position;
	}
}
void MultiResolutionHierarchy::construct_Es(const MatrixXu &mE, std::vector<tuple_E> &mEs){
	mEs.resize(mE.cols());

	for(int i=0;i<mE.cols();i++){
		mEs[i] = std::make_tuple(mE(0,i), mE(1,i), true, 0, Edge_tag::B, i, -1, 0);
	}
}
void MultiResolutionHierarchy::composit_edges_colors(MatrixXf &Result_Vs, std::vector<tuple_E> &Es_to_render, MatrixXf &Result_edges)
{
	Result_edges.resize(6, 2 * Es_to_render.size());
	//for rendering edges
	for (uint32_t i = 0; i < Es_to_render.size(); ++i) {
		Vector3f color;
		if (std::get<4>(Es_to_render[i]) == Edge_tag::R)
			color = Vector3f(1, 0, 0);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::B)
			color = Vector3f(0, 0, 1);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::D)
			color = Vector3f(0, 1, 0);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::H)
		color = Vector3f(1, 1, 1);

		uint32_t i0 = std::get<0>(Es_to_render[i]), i1 = std::get<1>(Es_to_render[i]);

		Result_edges.col(i * 2 + 0) << Result_Vs.col(i0), color;
		Result_edges.col(i * 2 + 1) << Result_Vs.col(i1), color;
	}
}
void MultiResolutionHierarchy::composit_edges_colors(Graph &g_, MatrixXf &Result_edges)
{
	Result_edges.resize(6, 2 * g_.Es.size());
	//for rendering edges
	for (uint32_t i = 0; i < g_.Es.size(); ++i)
	{
		auto &e = g_.Es[i];
		Vector3f color;
		if (e.color == Edge_tag::R)
			color = Vector3f(1, 0, 0);
		else if (e.color == Edge_tag::B)
			color = Vector3f(0, 0, 1);
		else if (e.color == Edge_tag::D)
			color = Vector3f(0, 1, 0);
		else if (e.color == Edge_tag::H)
			color = Vector3f(1, 1, 1);
		if(e.boundary)
			color = Vector3f(1, 1, 1);
		if (e.color != Edge_tag::B && (g_.Vs[e.vs[0]].singularity || g_.Vs[e.vs[1]].singularity))
			;// color = Vector3f(1, 1, 1);

		uint32_t i0 = e.vs[0], i1 = e.vs[1];

		// Result_edges.col(i * 2 + 0) << g_.Vs[i0].position, color;
		// Result_edges.col(i * 2 + 1) << g_.Vs[i1].position, color;
		Result_edges.col(i * 2 + 0) << g_.Vs[i0].v, color;
		Result_edges.col(i * 2 + 1) << g_.Vs[i1].v, color;
	}
}
void MultiResolutionHierarchy::edge_tagging2D() {

	composit_edges_colors(mV[0], mEs, E_rend);

	std::vector<std::vector<uint32_t>> PV_npes(mV[0].cols());
	for (uint32_t i = 0; i < mEs.size(); i++) {
		uint32_t v0 = std::get<0>(mEs[i]), v1 = std::get<1>(mEs[i]);
		PV_npes[v0].push_back(i); PV_npes[v1].push_back(i);
	}
	std::vector<bool> V_flag_(mV[0].cols(), false);
	std::vector<Vector3f> mO_center_;
	mV_flag.resize(mV[0].cols());
	mO_center.resize(3, mV[0].cols());
	while (true) {
		std::vector<uint32_t> v_pool, v_set;
		for (uint32_t i = 0; i < V_flag_.size(); i++) if (!V_flag_[i]) { v_pool.push_back(i); break; }
		if (!v_pool.size()) break;
		v_set = v_pool;
		V_flag_[v_pool[0]] = true;
		while (v_pool.size()) {
			std::vector<uint32_t> v_pool_sudo;
			for (uint32_t j = 0; j < v_pool.size(); j++)
				for (uint32_t k = 0; k < PV_npes[v_pool[j]].size(); k++) {
					uint32_t eid = PV_npes[v_pool[j]][k];
					uint32_t v0 = std::get<0>(mEs[eid]), v1 = std::get<1>(mEs[eid]);
					if (std::get<4>(mEs[eid]) == Edge_tag::R) {
						if (!V_flag_[v0]) v_pool_sudo.push_back(v0);
						if (!V_flag_[v1]) v_pool_sudo.push_back(v1);
					}
				}
			v_pool.clear();
			if (v_pool_sudo.size()) {
				v_pool.clear();
				for (uint32_t j = 0; j < v_pool_sudo.size(); j++) if (!V_flag_[v_pool_sudo[j]]) {
					v_pool.push_back(v_pool_sudo[j]); V_flag_[v_pool_sudo[j]] = true;
				}
				v_set.insert(v_set.end(), v_pool.begin(), v_pool.end());
			}
		}
		Vector3f center; center.setZero();
		std::vector<int> hard_constraints, boundary_vs;
		for (uint32_t j = 0; j < v_set.size(); j++)
		{
			if(BC[v_set[j]]) hard_constraints.push_back(v_set[j]);
			if(BV[0][v_set[j]]) boundary_vs.push_back(v_set[j]);
		}
		hard_constraints.clear();

		for (uint32_t j = 0; j < v_set.size(); j++){
				mV_flag[v_set[j]] = mO_center_.size();
		}

		if(hard_constraints.size())
		{
			center = mO[0].col(hard_constraints[0]);
		}else if(boundary_vs.size())
		{
			for (uint32_t j = 0; j < boundary_vs.size(); j++){
				center += mO[0].col(boundary_vs[j]);
			}
			center /= boundary_vs.size();			
		}else{
			for (uint32_t j = 0; j < v_set.size(); j++){
				center += mO[0].col(v_set[j]);
			}
			center /= v_set.size();
		}
		mO_center_.push_back(center);
		for (uint32_t j = 0; j < v_set.size(); j++)
			mO_center.col(v_set[j]) = center;
	}

	E_O_rend.resize(6, mEs.size() * 2);
	composit_edges_colors(mO_center, mEs, E_O_rend);

	E_I_rend = E_rend;
//mV_final
	mV_final.resize(3, mO_center_.size());
	for(int i=0;i<mO_center_.size();i++)
		mV_final.col(i) = mO_center_[i];
	
	std::vector<std::vector<int>> Es;
	std::vector<int> an_e(2);
	for(auto &e:mEs)if(std::get<4>(e) == Edge_tag::D)continue;
	else {
		an_e[0]=mV_flag[std::get<0>(e)];
		an_e[1]=mV_flag[std::get<1>(e)];
		if(an_e[0]>an_e[1])std::swap(an_e[0], an_e[1]);
		if(an_e[0]==an_e[1])continue;
		Es.push_back(an_e);
	}
	std::sort(Es.begin(),Es.end());
	Es.erase(std::unique(Es.begin(), Es.end()), Es.end());
	mE_final.resize(2,Es.size());
	for(int i=0;i<Es.size();i++){
		mE_final(0, i)=Es[i][0];
		mE_final(1, i)=Es[i][1];
	}
}
void MultiResolutionHierarchy::edge_tagging3D(){
	for (uint32_t i = 0; i < mEs.size(); ++i) {
		uint32_t v0 = std::get<0>(mEs[i]), v1 = std::get<1>(mEs[i]);
		//Quaternion q_next = Quaternion::applyRotation(mQ[0].col(v1), mQ[0].col(v0));	
		const Vector3f scale0 = mS[0].col(v0).head(3) * mScale;
        //Vector2f scale0 (S(2, j) * mScale, S(2, j) * mScale);
        //const Vector3f Invscale0 (1.0/scale0[0], 1.0/scale0[1], 1.0/scale0[2]);
		const Vector3f scale1 = mS[0].col(v1).head(3) * mScale;
        //Vector2f scale1 (S(2, j) * mScale, S(2, j) * mScale);
        //const Vector3f Invscale1 (1.0/scale1[0], 1.0/scale1[1], 1.0/scale1[2]);

		 auto a_pair = posy3D_completeInfo(mO[0].col(v0), mQ[0].col(v0),scale0, mO[0].col(v1), mQ[0].col(v1),scale1);
		 std::get<4>(mEs[i]) = get<0>(a_pair);
		 std::get<3>(mEs[i]) = get<1>(a_pair);
	}
	// for (uint32_t i = 0; i < mpFes.size(); ++i) {
	// 	short n_R = 0;
	// 	for (uint32_t j = 0; j < 3; j++) if (std::get<4>(mEs[mpFes[i][j]]) == Edge_tag::R) n_R++;
	// 	if (n_R == 2)
	// 		for (uint32_t j = 0; j < 3; j++) std::get<4>(mEs[mpFes[i][j]]) = Edge_tag::R;
	// }

	composit_edges_colors(mV[0], mEs, E_rend);

	std::vector<std::vector<uint32_t>> PV_npes(mV[0].cols());
	for (uint32_t i = 0; i < mEs.size(); i++) {
		uint32_t v0 = std::get<0>(mEs[i]), v1 = std::get<1>(mEs[i]);
		PV_npes[v0].push_back(i); PV_npes[v1].push_back(i);
	}
	std::vector<bool> V_flag_(mV[0].cols(),false);
	vector<vector<uint32_t>> VSets;
	std::vector<Vector3f> mO_center_;
	mV_flag.resize(mV[0].cols());
	mO_center.resize(3, mV[0].cols());
	while (true) {
		std::vector<uint32_t> v_pool, v_set;
		for (uint32_t i = 0; i < V_flag_.size(); i++) if (!V_flag_[i]) { v_pool.push_back(i); break; }
		if (!v_pool.size()) break;
		v_set = v_pool;
		V_flag_[v_pool[0]] = true;
		while (v_pool.size()) {
			std::vector<uint32_t> v_pool_sudo;
			for (uint32_t j = 0; j < v_pool.size(); j++)
				for (uint32_t k = 0; k < PV_npes[v_pool[j]].size(); k++) {
					uint32_t eid = PV_npes[v_pool[j]][k];
					uint32_t v0 = std::get<0>(mEs[eid]), v1 = std::get<1>(mEs[eid]);
					if (std::get<4>(mEs[eid]) == Edge_tag::R) {
						if (!V_flag_[v0]) v_pool_sudo.push_back(v0);
						if (!V_flag_[v1]) v_pool_sudo.push_back(v1);
					}
				}
			v_pool.clear();
			if (v_pool_sudo.size()) {
				for (uint32_t j = 0; j < v_pool_sudo.size(); j++) if (!V_flag_[v_pool_sudo[j]]) {
					v_pool.push_back(v_pool_sudo[j]); V_flag_[v_pool_sudo[j]] = true;
				}
				v_set.insert(v_set.end(), v_pool.begin(), v_pool.end());
			}
		}
		//cout<<"v_set: "<<v_set.size()<<endl;
		Vector3f center; center.setZero();
		for (uint32_t j = 0; j < v_set.size(); j++){
			mV_flag[v_set[j]] = mO_center_.size();
			center += mO[0].col(v_set[j]);
		}
		center /= v_set.size();
		mO_center_.push_back(center);
		for (uint32_t j = 0; j < v_set.size(); j++)
			mO_center.col(v_set[j]) = center;
		VSets.push_back(v_set);
	}

	E_O_rend.resize(6, mEs.size() * 2);
	composit_edges_colors(mO_center, mEs, E_O_rend);

	E_I_rend = E_rend;

//mV_final
	mV_final.resize(3, mO_center_.size());
	for(int i=0;i<mO_center_.size();i++)
		mV_final.col(i) = mO_center_[i];
	
	std::vector<std::vector<int>> Es;
	std::vector<int> an_e(2);
	for(auto &e:mEs)if(std::get<4>(e) == Edge_tag::D  || std::get<4>(e)==Edge_tag::H)continue;
	else {
		an_e[0]=mV_flag[std::get<0>(e)];
		an_e[1]=mV_flag[std::get<1>(e)];
		if(an_e[0]>an_e[1])std::swap(an_e[0], an_e[1]);
		if(an_e[0]==an_e[1])continue;
		Es.push_back(an_e);
	}
	std::sort(Es.begin(),Es.end());
	Es.erase(std::unique(Es.begin(), Es.end()), Es.end());
	mE_final.resize(2,Es.size());
	for(int i=0;i<Es.size();i++){
		mE_final(0, i)=Es[i][0];
		mE_final(1, i)=Es[i][1];
	}

}

void MultiResolutionHierarchy::edge_tagging(Graph &g_) {
	
	int new_color_r=0, new_color_b=0, new_color_d=0, new_color_h=0;
	for (auto &e : g_.Es) {
		int v0 = e.vs[0], v1 = e.vs[1];

		auto a_tuple = posy_info(g_, v0, v1);

		//if (step != 0 && std::get<0>(a_tuple) == Edge_tag::D)
		//{
		//	auto vv = get<2>(a_tuple);
		//	auto x = vv[0], y = vv[1];// , z = vv[2];
		//	if (x > y) swap(x, y);
		//	//if (y > z) swap(y, z);
		//	//if (x > y) swap(x, y);
		//	if (x / y < 0.5) {
		//		std::get<0>(a_tuple) = Edge_tag::B;
		//		cout << vv << endl;
		//		cout << "v0 v1: " << v0 << " " << v1 <<" "<<x/y<< endl;
		//		new_color_b++;
		//	}
		//}

		 //if(step!=0 && e.color != std::get<0>(a_tuple))
		 //{
		 //	switch (get<0>(a_tuple)) 
		 //			{
		 //				case Edge_tag::R: new_color_r++; 
		 //				if((g_.Vs[v0].position-g_.Vs[v1].position).norm() > 2*mScale)
		 //				{
		 //					Quaternion d0 =g_.Vs[v0].direction, d1 =g_.Vs[v1].direction;
		 //					cout<<"dirction: "<<d0.toMatrix()<<" "<<d1.toMatrix()<<endl;
		 //					cout<<"position: "<<g_.Vs[v0].position<<" "<<g_.Vs[v1].position<<endl;
		 //					cout<<"scale: "<<g_.Vs[v0].scale<<" "<<g_.Vs[v1].scale<<endl;
		 //				}
		 //				break;
		 //				case Edge_tag::B: new_color_b++; break;
		 //				case Edge_tag::D: new_color_d++; break;
		 //				case Edge_tag::H: new_color_h++; break;
		 //			}
		 //}

		e.color = std::get<0>(a_tuple);
		e.w = std::get<1>(a_tuple);
		e.posy_dif = std::get<2>(a_tuple);
	}
	cout<<new_color_r<<" "<<new_color_b<<" "<<new_color_d<<" "<<new_color_h<<" "<<endl;
	cout << mScale<<endl;
	step++;
}
std::tuple<Edge_tag, Float, VectorXf> MultiResolutionHierarchy::posy_info(Graph &g_, int v0, int v1, bool averaging)
{
	if(g_.dimension == 2)
	{
		Vector3f q_cur = g_.Vs[v0].direction;
		Vector3f q_next = g_.Vs[v1].direction;

		Vector2f scale0 = g_.Vs[v0].scale.head(2) * mScale;
		Vector2f scale1 = g_.Vs[v1].scale.head(2) * mScale;

		auto p0 = g_.Vs[v0].position;
		auto p1 = g_.Vs[v1].position;

		auto n0 = g_.Vs[v0].normal;
		auto n1 = g_.Vs[v1].normal;

		auto c0 = g_.Vs[v0].v;
		auto c1 = g_.Vs[v1].v;
	
		if(averaging)
		{
			double r0 = 1.0/g_.Reverse_V_map[v0].size();
			double r1 = 1.0/g_.Reverse_V_map[v1].size();

			q_cur.normalize();
			q_next.normalize();

			scale0 *= r0;
			scale1 *= r1;

			p0 *= r0;
			p1 *= r1;

			n0.normalize();
			n1.normalize();

			c0 *= r0;
			c1 *= r1;
		}

		return posy2D_completeInfo(p0, q_cur, n0, c0, scale0,
		p1, q_next, n1, c1, scale1, is_anisotropy);
	}
	else if(g_.dimension == 3)
	{
		Quaternion q_cur = g_.Vs[v0].direction;
		Quaternion q_next = g_.Vs[v1].direction;

		Vector3f scale0 = g_.Vs[v0].scale * mScale;
		Vector3f scale1 = g_.Vs[v1].scale * mScale;

		auto p0 = g_.Vs[v0].position;
		auto p1 = g_.Vs[v1].position;

		if(averaging)
		{
			double r0 = 1.0/g_.Reverse_V_map[v0].size();
			double r1 = 1.0/g_.Reverse_V_map[v1].size();

			q_cur.normalize();
			q_next.normalize();

			scale0 *= r0;
			scale1 *= r1;

			p0 *= r0;
			p1 *= r1;
		}

		return posy3D_completeInfo(p0, q_cur, scale0, p1, q_next, scale1, is_anisotropy);
	}
}
std::tuple<Edge_tag, Float, VectorXf> MultiResolutionHierarchy::posy_info(int dim, Graph_Node &v0, Graph_Node &v1, int vn0, int vn1, bool averaging)
{
	if (dim == 2)
	{
		Vector3f q_cur = v0.direction;
		Vector3f q_next = v1.direction;

		Vector2f scale0 = v0.scale.head(2) * mScale;
		Vector2f scale1 = v1.scale.head(2) * mScale;

		auto p0 = v0.position;
		auto p1 = v1.position;

		auto n0 = v0.normal;
		auto n1 = v1.normal;

		auto c0 = v0.v;
		auto c1 = v1.v;

		if (averaging)
		{
			double r0 = 1.0 / vn0;
			double r1 = 1.0 / vn1;

			q_cur.normalize();
			q_next.normalize();

			scale0 *= r0;
			scale1 *= r1;

			p0 *= r0;
			p1 *= r1;

			n0.normalize();
			n1.normalize();

			c0 *= r0;
			c1 *= r1;
		}

		return posy2D_completeInfo(p0, q_cur, n0, c0, scale0,
			p1, q_next, n1, c1, scale1, is_anisotropy);
	}
	else if (dim== 3)
	{
		Quaternion q_cur = v0.direction;
		Quaternion q_next = v1.direction;

		Vector3f scale0 = v0.scale * mScale;
		Vector3f scale1 = v1.scale * mScale;

		auto p0 = v0.position;
		auto p1 = v1.position;

		if (averaging)
		{
			double r0 = 1.0 / vn0;
			double r1 = 1.0 / vn1;

			q_cur.normalize();
			q_next.normalize();

			scale0 *= r0;
			scale1 *= r1;

			p0 *= r0;
			p1 *= r1;
		}

		return posy3D_completeInfo(p0, q_cur, scale0, p1, q_next, scale1, is_anisotropy);
	}
}
std::tuple<Edge_tag, Float, VectorXf> MultiResolutionHierarchy::posy_info(Graph &g_, int v0, int v1, Vector3f p0, Vector3f p1)
{
	if(g_.dimension == 2)
	{
		Vector3f q_cur = g_.Vs[v0].direction;
		Vector3f q_next = g_.Vs[v1].direction;

		const Vector2f scale0 = g_.Vs[v0].scale.head(2) * mScale;
		const Vector2f scale1 = g_.Vs[v1].scale.head(2) * mScale;
		return posy2D_completeInfo(p0, q_cur, g_.Vs[v0].normal, g_.Vs[v0].v, scale0,
		 p1, q_next, g_.Vs[v1].normal, g_.Vs[v1].v, scale1, is_anisotropy);
	}
	else if(g_.dimension == 3)
	{
		const Quaternion q_cur = g_.Vs[v0].direction;
		const Quaternion q_next = g_.Vs[v1].direction;

		const Vector3f scale0 = g_.Vs[v0].scale.head(3) * mScale;
		const Vector3f scale1 = g_.Vs[v1].scale.head(3) * mScale;
		return posy3D_completeInfo(p0, q_cur, scale0, p1, q_next, scale1, is_anisotropy);
	}
}
void MultiResolutionHierarchy::interpolate_twovs(VectorXf q0, Vector3f g0, Vector3f n0, Vector3f v0, Vector3f s0,
	VectorXf qj, Vector3f gj, Vector3f nj, Vector3f vj, Vector3f sj,
	VectorXf &qn, Vector3f &gn, Vector3f &nn, Vector3f &vn, Vector3f &sn,
	int dimension)
{
	gn = (g0 + gj) * 0.5;
	nn = (n0 + nj).normalized();
	vn = (v0 + vj) * 0.5;

	if(dimension == 2)
	{
		qn = (q0 + applyRotation(q0, n0, qj, nj)).normalized();
		
		if(is_anisotropy)
		{
			int which = findRotation(q0, n0, qj, nj);
			if (which == 1 || which == 3)
				std::swap(sj[0], sj[1]);
		}
		sn = (s0 + sj) * 0.5;		
	}else
	{
		Quaternion q1_ =Quaternion::applyRotation(qj, q0);
		qn = (q0 + q1_).normalized();

		Vector3f scale = sj;
		if(is_anisotropy)
		{
			Quaternion q = qj;
			Eigen::Matrix<Float, 3, 3> M0, M1;
			M0 = q.toMatrix();
			M1 = q1_.toMatrix();
			for(int j=0;j<3;j++){
				int id=0;Float dmax=-1;
				for(int k=0;k<3;k++){
					Float d = std::abs(M0.col(j).dot(M1.col(k)));
					if(d>dmax){dmax = d; id = k;}	
				}
				scale[j] = sj[id];
			}
		}
		sn = (s0 + scale) * 0.5;		
	}
}
void MultiResolutionHierarchy::detect_parallel_redundant_edges(Graph &g_, vector<int> &ledges)
{
	ledges.clear();
	double threshold = 1.5;

	for (const auto &v : g_.Vs)
	{
		auto &nes = v.nes;
		vector<int> es;
		for (auto &eid : nes)
			es.push_back(eid);

		for (int j = 0; j < es.size(); j++) 
		{
			for (int k = j + 1; k < es.size(); k++) 
			{
				if (g_.Es[es[j]].color != Edge_tag::B || g_.Es[es[k]].color != Edge_tag::B)
					continue;
				int v0 = g_.Es[es[j]].vs[0];
				int v1 = g_.Es[es[k]].vs[0];
				if (v0 == v.id)
					v0 = g_.Es[es[j]].vs[1];
				if (v1 == v.id)
					v1 = g_.Es[es[k]].vs[1];

				auto p0 = posy_info(g_, v.id, v0);
				auto p1 = posy_info(g_, v.id, v1);

				double len0, len1;
				if (get<2>(p0)[0] > 0.5 && get<2>(p1)[0] > 0.5)
				{
					len0 = get<2>(p0)[0];
					len1 = get<2>(p1)[0];
				}
				else if (get<2>(p0)[1] > 0.5 && get<2>(p1)[1] > 0.5)
				{
					len0 = get<2>(p0)[1];
					len1 = get<2>(p1)[1];
				}
				else
					continue;

				auto &q0 = g_.Vs[v.id].direction;
				auto &g0 = g_.Vs[v.id].position;
				auto &n0 = g_.Vs[v.id].normal;
				auto &v0_ = g_.Vs[v.id].v;
				VectorXf s0 = g_.Vs[v.id].scale;
				auto &qj = g_.Vs[v0].direction;
				auto &gj = g_.Vs[v0].position;
				auto &nj = g_.Vs[v0].normal;
				auto &vj = g_.Vs[v0].v;
				Vector3f sj = g_.Vs[v0].scale;
				auto &qk = g_.Vs[v1].direction;
				auto &gk = g_.Vs[v1].position;
				auto &nk = g_.Vs[v1].normal;
				auto &vk = g_.Vs[v1].v;
				Vector3f sk = g_.Vs[v1].scale;

				VectorXf qn;
				Vector3f gn;
				Vector3f nn;
				Vector3f vn;
				Vector3f sn;
				VectorXf isn;
				std::tuple<Edge_tag, Float, VectorXf> a_posy;

				if (std::round(len0 / len1) >= threshold)
				{
					interpolate_twovs(q0, g0, n0, v0_, s0, qj, gj, nj, vj, sj, qn, gn, nn, vn, sn, g_.dimension);

					std::tuple<Edge_tag, Float, VectorXf> a_posy;
					if(g_.dimension == 2)
					{
						a_posy = posy2D_completeInfo(
							gk, qk, nk, vk, sk.head(2)*mScale,
							gn, qn, nn, vn, sn.head(2)*mScale, is_anisotropy
						);
					}else
					{
						a_posy = posy3D_completeInfo(
							gk, qk, sk*mScale, 
							gn, qn, sn*mScale, is_anisotropy);
					}
				
					if (std::get<0>(a_posy) != Edge_tag::R) continue;

					ledges.push_back(es[j]);
					ledges.push_back(v1);
				}
				else if (std::round(len1 / len0) >= threshold)
				{
					interpolate_twovs(q0, g0, n0, v0_, s0, qk, gk, nk, vk, sk, qn, gn, nn, vn, sn, g_.dimension);

					if(g_.dimension == 2)
					{
						a_posy = posy2D_completeInfo(
							gj, qj, nj, vj, sj.head(2)*mScale,
							gn, qn, nn, vn, sn.head(2)*mScale, is_anisotropy
						);
					}else
					{
						a_posy = posy3D_completeInfo(
							gj, qj, sj*mScale,
							gn, qn, sn*mScale, is_anisotropy);
					}

					if (std::get<0>(a_posy) != Edge_tag::R) continue;

					ledges.push_back(es[k]);
					ledges.push_back(v0);
				}
			}
		}
	}
}
void MultiResolutionHierarchy::split_parallel_redundant_edges(Graph &g_, vector<int> &ledges)
{
	std::vector<bool> E_tag(g_.Es.size(), true);
	for (int i = 0; i < ledges.size()/2; i++) {
		auto & e = g_.Es[ledges[2*i]];

		Graph_Node gn;
		gn.id = g_.Vs.size();
		gn.boundary = e.boundary;

		auto q0 = g_.Vs[e.vs[0]].direction;
		auto q1 = g_.Vs[e.vs[1]].direction;
		auto s0 = g_.Vs[e.vs[0]].scale;
		auto s1 = g_.Vs[e.vs[1]].scale;
		
		interpolate_twovs(q0,g_.Vs[e.vs[0]].position,g_.Vs[e.vs[0]].normal,g_.Vs[e.vs[0]].v,s0, 
							q1, g_.Vs[e.vs[1]].position,g_.Vs[e.vs[1]].normal,g_.Vs[e.vs[1]].v,s1,
							gn.direction,gn.position,gn.normal, gn.v, gn.scale, g_.dimension);

		g_.Vs.push_back(gn);

		std::vector<int> sharedvs;
		sharedvs.push_back(ledges[2 * i + 1]);
		sharedvs.push_back(e.vs[0]);
		sharedvs.push_back(e.vs[1]);
		for (auto &vid : sharedvs)
		{
			Graph_Edge ge;
			ge.vs.push_back(gn.id);
			ge.vs.push_back(vid);
			ge.boundary = false;
			if (vid == e.vs[0] || vid == e.vs[1])
				ge.boundary = e.boundary;
			ge.id = g_.Es.size();
			ge.time = 0;
			g_.Es.push_back(ge);
			E_tag.push_back(true);

			auto a_tuple = posy_info(g_, gn.id, vid);
			
			g_.Es[ge.id].color = std::get<0>(a_tuple);
			g_.Es[ge.id].w = std::get<1>(a_tuple);
			g_.Es[ge.id].posy_dif = std::get<2>(a_tuple);

			g_.Vs[gn.id].nvs.insert(vid);
			g_.Vs[vid].nvs.insert(gn.id);
			g_.Vs[gn.id].nes.insert(ge.id);
			g_.Vs[vid].nes.insert(ge.id);
		}
		g_.Vs[e.vs[0]].nvs.erase(e.vs[1]);
		g_.Vs[e.vs[1]].nvs.erase(e.vs[0]);
		g_.Vs[e.vs[0]].nes.erase(e.id);
		g_.Vs[e.vs[1]].nes.erase(e.id);

		E_tag[e.id] = false;
	}
	if (ledges.size()) {
		auto GEs = g_.Es;
		g_.Es.clear();
		std::vector<int> E_map(GEs.size(), -1);
		for (auto &e : GEs)
		{
			if (E_tag[e.id])
			{
				E_map[e.id] = g_.Es.size();
				e.id = g_.Es.size();
				g_.Es.push_back(e);
			}
		}
		for (auto &v : g_.Vs)
		{
			auto nes = v.nes;
			v.nes.clear();
			for (const auto &eid : nes)
				if (E_map[eid] != -1)v.nes.insert(E_map[eid]);
		}
	}
}
void MultiResolutionHierarchy::connect_close_points(Graph &g_)
{
	MatrixXf V;
	MatrixXu F;
	construct_G2VE(g_, V);

	BVH *bvh_local = new BVH(&F, &V, mAABB);
	bvh_local->build();

	for (auto &v:g_.Vs) 
	{
		double max_s = g_.Vs[v.id].scale.maxCoeff();
		//double radius = 0.5 * max_s * mScale;
		double radius = 0.1 * max_s * mScale;
		vector<uint32_t> result;
		bvh_local->findNearestWithRadius(v.position, radius, result);

		std::vector<int> newvs;
		for (auto &vid : result)
		{
			if (!v.nvs.count(vid))
				newvs.push_back(vid);
		}
		//cout<<"connected "<<newvs.size()<<endl;
		for (auto &vid : newvs)
		{
			Graph_Edge ge;
			ge.vs.push_back(v.id);
			ge.vs.push_back(vid);
			ge.boundary = false;
			ge.id = g_.Es.size();
			ge.time = 0;
			g_.Es.push_back(ge);

			auto a_tuple = posy_info(g_, v.id, vid);
			g_.Es[ge.id].color = std::get<0>(a_tuple);
			g_.Es[ge.id].w = std::get<1>(a_tuple);
			g_.Es[ge.id].posy_dif = std::get<2>(a_tuple);

			g_.Vs[v.id].nvs.insert(vid);
			g_.Vs[vid].nvs.insert(v.id);
			g_.Vs[v.id].nes.insert(ge.id);
			g_.Vs[vid].nes.insert(ge.id);
		}
	}

	delete bvh_local;
}
bool MultiResolutionHierarchy::collapse_fuse_edges(Graph &g_)
{
	struct LessThanGraph { bool operator()(const Graph_Edge& lhs, const Graph_Edge& rhs) const { return (lhs.w > rhs.w); } };

	std::priority_queue<Graph_Edge, std::vector<Graph_Edge>, LessThanGraph> GEs;

	for (auto & e:g_.Es)
		GEs.push(e);

	Graph G_ = g_;

	composit_edges_colors(g_, E_rend_o);

	std::vector<int> E_TimeStamp(g_.Es.size(), 0);
	std::vector<bool> V_Tag(g_.Vs.size(), true), E_Tag(g_.Es.size(), true);

	g_.V_map.clear();
	g_.V_map.resize(g_.Vs.size());
	for (int i = 0; i < g_.Vs.size(); i++)
		g_.V_map[i] = i;
	g_.Reverse_V_map.clear();
	g_.Reverse_V_map.resize(g_.Vs.size());
	for (int i = 0; i < g_.Vs.size(); i++)
		g_.Reverse_V_map[i].insert(i);

	g_.Reverse_E_map.clear();
	g_.Reverse_E_map.resize(g_.Es.size());
	for (int i = 0; i < g_.Es.size(); i++)
		g_.Reverse_E_map[i].insert(i);

	int count_passed = 0;
	while(!GEs.empty())
	{
		auto e = GEs.top();
		GEs.pop();
		e = g_.Es[e.id];
		if(!E_Tag[e.id]) 
			continue;
		if(e.time < E_TimeStamp[e.id]) 
			continue;
		
		int v0 = e.vs[0], v1 = e.vs[1];

		if(e.color == Edge_tag::B) continue;
		if(e.color == Edge_tag::R) 
		{//collapse
			if (g_.dimension == 2 && g_.Vs[v0].boundary && g_.Vs[v1].boundary && !e.boundary)
			{
				count_passed++;
				continue;//internal edge with boundary vs
			}
			
			if((!g_.Vs[v0].boundary && !g_.Vs[v1].boundary)|| (g_.Vs[v0].boundary && g_.Vs[v1].boundary))
			{
				bool is_boundary = (g_.Vs[v0].boundary && g_.Vs[v1].boundary);

				auto rvs = g_.Reverse_V_map[v0];
				rvs.insert(g_.Reverse_V_map[v1].begin(), g_.Reverse_V_map[v1].end());
				
				auto gn = g_.Vs[v1];
				gn.v = gn.position = gn.normal = Vector3f::Zero();
				gn.scale.setConstant(0);

				if(g_.dimension ==2)
					gn.direction = Vector3f::Zero();
				else
					gn.direction = Vector4f::Zero();

				auto q0 = gn.direction;
				auto n0 = gn.normal;
				auto s0 = gn.scale;
				for (auto &vid : rvs)
				{
					q0 = G_.Vs[vid].direction;
					n0 = G_.Vs[vid].normal;
					s0 = G_.Vs[vid].scale;
					break;
				}
				int nv = 0;
				for (auto &vid : rvs)
				{
					if (is_boundary)
					{
						if (!G_.Vs[vid].boundary)							
							continue;
					}
					nv++;
					gn.v += G_.Vs[vid].v;
					gn.position += G_.Vs[vid].position;
					gn.normal += G_.Vs[vid].normal;					
					auto s_1 = G_.Vs[vid].scale, s_1_ = s_1;

					if(g_.dimension == 2)
					{
						gn.direction += applyRotation(q0, n0, G_.Vs[vid].direction, G_.Vs[vid].normal);
						
						if(is_anisotropy)
						{
							int which = findRotation(q0, n0, G_.Vs[vid].direction, G_.Vs[vid].normal);
							if (which == 1 || which == 3)
								std::swap(s_1[0], s_1[1]);
						}
					}else
					{
						Quaternion q1_ =Quaternion::applyRotation(G_.Vs[vid].direction, q0);
						gn.direction += q1_;

						if(is_anisotropy)
						{
							Quaternion q = G_.Vs[vid].direction;
							Eigen::Matrix<Float, 3, 3> M0, M1;
							M0 = q.toMatrix();
							M1 = q1_.toMatrix();
							for(int j=0;j<3;j++){
								int id=0;Float dmax=-1;
								for(int k=0;k<3;k++){
									Float d = std::abs(M0.col(j).dot(M1.col(k)));
									if(d>dmax){dmax = d; id = k;}	
								}
								s_1[j] = s_1_[id];
							}
						}
					}
					gn.scale += s_1;
				}
				gn.v /= nv;
				gn.position /= nv;
				gn.normal.normalize();
				gn.direction.normalize();
				gn.scale /= nv;

				g_.Vs[v1] = gn;
			}
			else if(!g_.Vs[e.vs[0]].boundary && g_.Vs[e.vs[1]].boundary)
			{
				;
			}
			else if(g_.Vs[e.vs[0]].boundary && !g_.Vs[e.vs[1]].boundary)
			{
				std::swap(v0, v1);
			}
			
			//new v
			V_Tag[v0] = false;
			E_Tag[e.id] = false;
			//update connectivity
			auto & v = g_.Vs[v1];
			auto & v_del = g_.Vs[v0];
			v.boundary = v.boundary || v_del.boundary;

			g_.Reverse_V_map[v.id].insert(g_.Reverse_V_map[v_del.id].begin(), g_.Reverse_V_map[v_del.id].end());
			g_.Reverse_V_map[v_del.id].clear();
			for (auto &vid : g_.Reverse_V_map[v.id])
				g_.V_map[vid] = v.id;
				
			const auto &nvs0 = v.nvs;
			const auto &nvs1 = v_del.nvs;
			std::vector<int> sharedvs, commones;
			unordered_set_intersection_own(nvs0, nvs1, sharedvs);
			for(auto &vid: sharedvs)
			{
				const auto &nes0 = g_.Vs[vid].nes;
				const auto &nes1 = v_del.nes;
				std::vector<int> sharedes0;
				unordered_set_intersection_own(nes0, nes1, sharedes0);
				if(!sharedes0.size())
					std::cout<<"error"<<std::endl;

				const auto &nes2 = v.nes;
				std::vector<int> sharedes1;
				unordered_set_intersection_own(nes0, nes2, sharedes1);
				if (!sharedes1.size())
					std::cout << "error" << std::endl;

				int e0 = sharedes0[0];
				int e1 = sharedes1[0];

				g_.Es[e0].boundary = g_.Es[e0].boundary || g_.Es[e1].boundary;
				E_Tag[e1] = false;
				commones.push_back(e1);

				g_.Reverse_E_map[e0].insert(g_.Reverse_E_map[e1].begin(), g_.Reverse_E_map[e1].end());
				g_.Reverse_E_map[e1].clear();
			}
			//v.nvs	
			v.nvs.insert(v_del.nvs.begin(),v_del.nvs.end());
			v.nvs.erase(v_del.id);
			v.nvs.erase(v.id);

			for(auto &vid: v_del.nvs)
			{
				if(vid == v.id)
					continue;
				g_.Vs[vid].nvs.erase(v_del.id);
				g_.Vs[vid].nvs.insert(v.id);
			}
			//v.nes 
			v.nes.insert(v_del.nes.begin(),v_del.nes.end());	
			v.nes.erase(e.id);
			for (auto &eid : commones) {
				auto &vs = g_.Es[eid].vs;
				g_.Vs[vs[0]].nes.erase(eid);
				g_.Vs[vs[1]].nes.erase(eid);
			}

			for(auto & eid: v_del.nes)
			{
				if(!E_Tag[eid]) continue;			
				if(g_.Es[eid].vs[0] == v_del.id) g_.Es[eid].vs[0] = v.id;
				else if(g_.Es[eid].vs[1] == v_del.id) g_.Es[eid].vs[1] = v.id;
				else
				{
					std::cout<<"error "<<std::endl;
				}
			}
			for(auto &eid: v.nes)
			{
				int votes[4] = {0 ,0 ,0, 0};
				auto &reverse_nes = g_.Reverse_E_map[eid];
				for (auto & reid : reverse_nes)
				{
					v0 = G_.Es[reid].vs[0];
					v1 = G_.Es[reid].vs[1];

					int map_v0 = g_.V_map[v0];
					int map_v1 = g_.V_map[v1];

					auto &g0 = g_.Vs[map_v0].position;
					auto &gj = g_.Vs[map_v1].position;

					auto a_posy = posy_info(G_, v0, v1, g0, gj);

					switch (get<0>(a_posy)) 
					{
						case Edge_tag::R: votes[0]++; break;
						case Edge_tag::B: votes[1]++; break;
						case Edge_tag::D: votes[2]++; break;
						case Edge_tag::H: votes[3]++; break;
					}
				}
				
				auto color = (Edge_tag)std::distance(votes, std::max_element(votes, votes + 4));;
				if(g_.Es[eid].color != color)
				{
					g_.Es[eid].color = color;
					E_TimeStamp[eid]++;
					g_.Es[eid].time = E_TimeStamp[eid];
					GEs.push(g_.Es[eid]);
				}				
			}
			//hard constraint alignment
		}
		else
		{//dissolve
			//continue;
			if(g_.dimension==2 && e.boundary) 
				continue;
			E_Tag[e.id] = false;			
			g_.Vs[e.vs[0]].nes.erase(e.id);
			g_.Vs[e.vs[1]].nes.erase(e.id);
			g_.Vs[e.vs[0]].nvs.erase(e.vs[1]);
			g_.Vs[e.vs[1]].nvs.erase(e.vs[0]);
		}	
	}
	for (int i=0;i<g_.Es.size();i++)
	{
		auto &es = g_.Reverse_E_map[i];
		if (!es.size())
			continue;
		for (auto &eid : es)
			g_.Es[i].boundary = g_.Es[i].boundary || g_.Es[eid].boundary;
	}
	//rendering
	for (auto &v : g_.Vs)
	{
		v.position = g_.Vs[g_.V_map[v.id]].position;
	}
	composit_edges_colors(g_, E_O_rend_o);

	//re-indexing
	Graph g_temp;
	g_temp.dimension = g_.dimension;

	g_temp.Vs.reserve(g_.Vs.size());
	g_temp.Es.reserve(g_.Es.size());
	Eigen::VectorXi v_map(V_Tag.size());
	v_map.setConstant(-1);
	int vn = 0;	
	for (int i = 0; i < V_Tag.size(); i++)
	{
		if(!V_Tag[i]) continue;

		v_map[i] = vn;
		auto v = g_.Vs[i];
		v.id = vn++;
		v.nvs.clear();
		v.nes.clear();
		g_temp.Vs.push_back(v);
		g_temp.Reverse_V_map.push_back(g_.Reverse_V_map[i]);
	}

	int en=0;
	int count = 0;
	for (int i = 0; i < E_Tag.size(); i++)
	{
		if (!E_Tag[i])
		{
			count++;
			continue;
		}

		auto e = g_.Es[i];
		e.id = en++;
		e.vs[0] = v_map[g_.Es[i].vs[0]];
		e.vs[1] = v_map[g_.Es[i].vs[1]];

		if(e.vs[0]>e.vs[1]) 
			std::swap(e.vs[0], e.vs[1]);
		g_temp.Es.push_back(e);

		g_temp.Vs[e.vs[0]].nvs.insert(e.vs[1]);
		g_temp.Vs[e.vs[1]].nvs.insert(e.vs[0]);

		g_temp.Vs[e.vs[0]].nes.insert(e.id);
		g_temp.Vs[e.vs[1]].nes.insert(e.id);

		g_temp.Reverse_E_map.push_back(g_.Reverse_E_map[i]);
	}
	g_ = g_temp;

	std::cout << "passed " << count_passed << " edges" << std::endl;
	std::cout<<"removed "<<count<<" edges"<<std::endl;
	return count == 0? false : true;
}
void MultiResolutionHierarchy::remove_dangling(Graph &g_)
{
	int removed_nv=0, added_ne =0, removed_ne=0;

	std::vector<bool> V_Tag(g_.Vs.size(), true);
	for(auto &v: g_.Vs)
	{
		if(v.boundary)
		{
			if(v.nvs.size()==1)
			{
				V_Tag[v.id] = false;
				removed_nv++;
			}
			if(v.nvs.size() == 2)
			{
				V_Tag[v.id] = false;
				removed_nv++;
				added_ne++;
			}
		} 
	}
	for(auto &v: g_.Vs)
	{
		if(v.boundary)
		{
			if(v.nvs.size() == 2)
			{
				Graph_Edge e;
				e.id = g_.Es.size();
				for(auto &vid: v.nvs)
					e.vs.push_back(vid);
				assert(e.vs.size()==2);

				g_.Es.push_back(e);
			}
		}
	}

	std::vector<bool> E_Tag(g_.Es.size(), true);
	for(auto &e: g_.Es)
	{
		E_Tag[e.id] = V_Tag[e.vs[0]]&&V_Tag[e.vs[1]];
		if(!E_Tag[e.id])
			removed_ne++;
	}

	//re-indexing
	Graph g_temp;
	g_temp.dimension = g_.dimension;

	g_temp.Vs.reserve(g_.Vs.size());
	g_temp.Es.reserve(g_.Es.size());
	Eigen::VectorXi v_map(V_Tag.size());
	v_map.setConstant(-1);
	int vn = 0;	
	for (int i = 0; i < V_Tag.size(); i++)
	{
		if(!V_Tag[i]) 
			continue;

		v_map[i] = vn;
		auto v = g_.Vs[i];
		v.id = vn++;
		v.nvs.clear();
		v.nes.clear();
		g_temp.Vs.push_back(v);
	}

	int en=0;
	for (int i = 0; i < E_Tag.size(); i++)
	{
		if (!E_Tag[i])
			continue;
		auto e = g_.Es[i];
		e.id = en++;
		e.vs[0] = v_map[g_.Es[i].vs[0]];
		e.vs[1] = v_map[g_.Es[i].vs[1]];

		if(e.vs[0]>e.vs[1]) 
			std::swap(e.vs[0], e.vs[1]);
		g_temp.Es.push_back(e);

	}
	g_ = g_temp;

	std::cout<<"removed v e"<<removed_nv<<" "<<removed_ne<<std::endl;	
	std::cout<<"added e"<<added_ne<<std::endl;	
}
void MultiResolutionHierarchy::graph_extraction()
{
	//construct_Graph(m, g);
	bool changes = true;
	int iter = 0;
	edge_tagging(g);
	int vn_pre = g.Vs.size(), en_pre = g.Es.size();
	while(changes)
	{
		std::cout<<iter++<<" iteration..."<<std::endl;
		changes = false;

		connect_close_points(g);

		vector<int> ledges;
		detect_parallel_redundant_edges(g, ledges);
		split_parallel_redundant_edges(g, ledges);
		cout <<"split "<< ledges.size()/2 << endl;
		//changes = split_long_edges(g);
		changes	= collapse_fuse_edges(g);

		//edge_tagging(g);
		if (vn_pre == g.Vs.size() && en_pre == g.Es.size())
			changes = false;
		else
		{
			vn_pre = g.Vs.size();
			en_pre = g.Es.size();
		}
	}
	//remove_dangling(g);
}

void MultiResolutionHierarchy::group_representative(const Graph &g_, const vector<int> &vs, Graph_Node &gn, bool averaging)
{
	if(!vs.size()) return;
	gn.v = gn.position = gn.normal = Vector3f::Zero();
	gn.scale.setConstant(0);

	if(g_.dimension ==2)
		gn.direction = Vector3f::Zero();
	else
		gn.direction = Vector4f::Zero();

	auto q0 = g_.Vs[vs[0]].direction;
	auto n0 = g_.Vs[vs[0]].normal;
	auto s0 = g_.Vs[vs[0]].scale;
	
	vector<int> bvs;
	for (auto &vid:vs) if(g_.Vs[vid].boundary) bvs.push_back(vid);
	
	int n = vs.size();
	gn.boundary = false;
	if(bvs.size())
	{
		gn.boundary = true;
		n = bvs.size();
	}

	for (int i =0;i<n;i++)
	{
		int vid = -1;
		if(gn.boundary) vid = bvs[i];
		else vid = vs[i];

		gn.v += g_.Vs[vid].v;
		gn.position += g_.Vs[vid].position;
		gn.normal += g_.Vs[vid].normal;					
		auto s_1 = g_.Vs[vid].scale, s_1_ = s_1;

		if(g_.dimension == 2) {
			gn.direction += applyRotation(q0, n0, g_.Vs[vid].direction, g_.Vs[vid].normal);
			
			if(is_anisotropy)
			{
				int which = findRotation(q0, n0, g_.Vs[vid].direction, g_.Vs[vid].normal);
				if (which == 1 || which == 3) std::swap(s_1[0], s_1[1]);
			}
		}else {
			Quaternion q1_ =Quaternion::applyRotation(g_.Vs[vid].direction, q0);
			gn.direction += q1_;

			if(is_anisotropy) {
				Quaternion q = g_.Vs[vid].direction;
				Eigen::Matrix<Float, 3, 3> M0, M1;
				M0 = q.toMatrix();
				M1 = q1_.toMatrix();
				for(int j=0;j<3;j++){
					int id=0;Float dmax=-1;
					for(int k=0;k<3;k++){
						Float d = std::abs(M0.col(j).dot(M1.col(k)));
						if(d>dmax){dmax = d; id = k;}	
					}
					s_1[j] = s_1_[id];
				}
			}
		}
		gn.scale += s_1;
	}
	if(averaging)
	{
		gn.v /= n;
		gn.position /= n;
		gn.normal.normalize();
		gn.direction.normalize();
		gn.scale /= n;
	}
}
void MultiResolutionHierarchy::split_group(Graph &g_, vector<vector<int>> & sets)
{
	struct LessThanGraph { bool operator()(const Graph_Edge& lhs, const Graph_Edge& rhs) const { return (lhs.w > rhs.w); } };

	std::priority_queue<Graph_Edge, std::vector<Graph_Edge>, LessThanGraph> GEs;

	for (auto & e:g_.Es) GEs.push(e);

	std::vector<int> E_TimeStamp(g_.Es.size(), 0);
	std::vector<bool> V_Tag(g_.Vs.size(), true), E_Tag(g_.Es.size(), true);

	g_.V_map.clear();
	g_.V_map.resize(g_.Vs.size());
	for (int i = 0; i < g_.Vs.size(); i++) g_.V_map[i] = i;
	g_.Reverse_V_map.clear();
	g_.Reverse_V_map.resize(g_.Vs.size());
	for (int i = 0; i < g_.Vs.size(); i++) g_.Reverse_V_map[i].insert(i);

	g_.Reverse_E_map.clear();
	g_.Reverse_E_map.resize(g_.Es.size());
	for (int i = 0; i < g_.Es.size(); i++) g_.Reverse_E_map[i].insert(i);

	vector<int> denometor(g_.Vs.size(), 1);
	while(!GEs.empty())
	{
		auto e = GEs.top();
		GEs.pop();
		e = g_.Es[e.id];
		if(!E_Tag[e.id]) 
			continue;
		if(e.time < E_TimeStamp[e.id]) 
			continue;
		
		int v0 = e.vs[0], v1 = e.vs[1];
		v0 = g_.V_map[v0];
		v1 = g_.V_map[v1];

		if(v0 == v1) continue;
		if(e.color != Edge_tag::R) continue;	
		E_Tag[e.id] = false;

		//update connectivity
		auto & v = g_.Vs[v1];
		auto & v_del = g_.Vs[v0];

		g_.Reverse_V_map[v.id].insert(g_.Reverse_V_map[v_del.id].begin(), g_.Reverse_V_map[v_del.id].end());
		g_.Reverse_V_map[v_del.id].clear();
		for (auto &vid : g_.Reverse_V_map[v.id])
			g_.V_map[vid] = v.id;
			
		const auto &nvs0 = v.nvs;
		const auto &nvs1 = v_del.nvs;
		std::vector<int> sharedvs, commones;
		unordered_set_intersection_own(nvs0, nvs1, sharedvs);
		for(auto &vid: sharedvs)
		{
			const auto &nes0 = g_.Vs[vid].nes;
			const auto &nes1 = v_del.nes;
			std::vector<int> sharedes0;
			unordered_set_intersection_own(nes0, nes1, sharedes0);
			if(!sharedes0.size())
				std::cout<<"error"<<std::endl;

			const auto &nes2 = v.nes;
			std::vector<int> sharedes1;
			unordered_set_intersection_own(nes0, nes2, sharedes1);
			if (!sharedes1.size())
				std::cout << "error" << std::endl;

			int e0 = sharedes0[0];
			int e1 = sharedes1[0];

			E_Tag[e1] = false;
			commones.push_back(e1);

			g_.Reverse_E_map[e0].insert(g_.Reverse_E_map[e1].begin(), g_.Reverse_E_map[e1].end());
			g_.Reverse_E_map[e1].clear();
		}
		//v.nvs	
		v.nvs.insert(v_del.nvs.begin(),v_del.nvs.end());
		v.nvs.erase(v_del.id);
		v.nvs.erase(v.id);

		for(auto &vid: v_del.nvs)
		{
			if(vid == v.id)
				continue;
			g_.Vs[vid].nvs.erase(v_del.id);
			g_.Vs[vid].nvs.insert(v.id);
		}
		//v.nes 
		v.nes.insert(v_del.nes.begin(),v_del.nes.end());	
		v.nes.erase(e.id);
		for (auto &eid : commones) {
			auto &vs = g_.Es[eid].vs;
			g_.Vs[vs[0]].nes.erase(eid);
			g_.Vs[vs[1]].nes.erase(eid);
		}

		Graph_Node gn;
		vector<int> vs;
		vs.push_back(v0);
		vs.push_back(v1);		
		group_representative(g_, vs, gn, false);
		g_.Vs[v1].position = gn.position;
		g_.Vs[v1].direction = gn.direction;
		g_.Vs[v1].scale = gn.scale;
		g_.Vs[v1].v = gn.v;
		g_.Vs[v1].normal = gn.normal;

		vector<int> bvs;
		for (auto &vid : g_.Reverse_V_map[v1]) if (g_.Vs[vid].boundary) bvs.push_back(vid);

		denometor[v1] = vs.size();
		gn.boundary = false;
		if (bvs.size()) denometor[v1] = vs.size();
		else denometor[v1] = bvs.size();

		Graph_Node gm;
		for(auto &eid: v.nes)
		{
			auto v0_ = g_.Es[eid].vs[0];
			auto v1_ = g_.Es[eid].vs[1];
			v0_ = g_.V_map[v0_];
			v1_ = g_.V_map[v1_];
			
			int vid = v0_;
			if (v0_ == v1) {
				vid = v1_;
				gm = g_.Vs[v1_];
			}
			else {
				gm = g_.Vs[v0_];
			}
			auto a_posy = posy_info(g_.dimension, gn, gm, denometor[v1], denometor[vid], true);
			auto color = get<0>(a_posy);
			if (color == Edge_tag::R && get<1>(a_posy) > 0.1)
				color = Edge_tag::B;
			if(g_.Es[eid].color != color)
			{
				g_.Es[eid].color = color;
				g_.Es[eid].w = get<1>(a_posy);
				E_TimeStamp[eid]++;
				g_.Es[eid].time = E_TimeStamp[eid];
				GEs.push(g_.Es[eid]);
			}
		}
	}
	//re-grouping
	sets.clear();
	vector<bool> v_flag(g_.Vs.size(),true);
	while (true) {		
		vector<int> v_set;
		int id = std::find(v_flag.begin(), v_flag.end(), true) - v_flag.begin();
		if(id >= v_flag.size())
			break;
		id = g_.V_map[id];
		for(int i=0;i<g_.V_map.size();i++){
			if(g_.V_map[i] == id){
				v_set.push_back(i);
				v_flag[i] = false;
			}
		}
		sets.push_back(v_set);
	}
}
void MultiResolutionHierarchy::split_group_based_on_color(Graph &g_, vector<vector<int>> & sets)
{
	for(auto & e: g_.Es)
	{
		if(e.color != Edge_tag::R)
		{
			std::vector<bool> v_flag(g_.Vs.size(),true);
			std::vector<bool> e_flag(g_.Es.size(),true);
			vector<vector<int>> vsets;

			auto &color = e.color;
			e.color = Edge_tag::R;

			e_flag[e.id] = false;
			while (true) {

				int id = std::find(v_flag.begin(), v_flag.end(), true) - v_flag.begin();
				if(id >= v_flag.size())
					break;

				std::vector<int> v_set;
				std::queue<int> pool;
				pool.push(id);
				v_flag[id] = true;

				while (!pool.empty()) {
					auto id_ = pool.front();
					v_set.push_back(id_);

					pool.pop();

					for (const auto &eid: g_.Vs[id_].nes) {
						if (e_flag[eid] && g_.Es[eid].color == Edge_tag::R) {
							e_flag[eid] = false;

							auto &vs = g_.Es[eid].vs;
							if (v_flag[vs[0]]) 
							{
								pool.push(vs[0]);
								v_flag[vs[0]] = false;
							}
							if (v_flag[vs[1]]) 
							{
								pool.push(vs[1]);
								v_flag[vs[1]] = false;
							}
						}
					}
				}
				vsets.push_back(v_set);
			}
			e.color = color;
			if(vsets.size()==2)
			{
				sets = vsets;
				return;
			}
		}
	}
}
void MultiResolutionHierarchy::grouping_collapse(Graph &g_)
{
	//grouping
	std::vector<int> v_flag(g_.Vs.size(),-1);
	std::vector<bool> e_flag(g_.Es.size(),true);
	vector<vector<int>> vsets;
	int ng = 0;
	while (true) {		
		std::vector<int> v_set;
		std::queue<int> pool;

		int id = std::find(v_flag.begin(), v_flag.end(), -1) - v_flag.begin();
		if(id >= v_flag.size())
			break;
		
		pool.push(id);
		v_flag[id] = ng;

		while (!pool.empty()) {
			auto id_ = pool.front();
			v_set.push_back(id_);

			pool.pop();

			for (const auto &eid: g_.Vs[id_].nes) {
				if (e_flag[eid] && g_.Es[eid].color == Edge_tag::R) {
					e_flag[eid] = false;

					auto &vs = g_.Es[eid].vs;
					if (v_flag[vs[0]] == -1) 
					{
						pool.push(vs[0]);
						v_flag[vs[0]] = ng;
					}
					if (v_flag[vs[1]] == -1) 
					{
						pool.push(vs[1]);
						v_flag[vs[1]] = ng;
					}
				}
			}
		}
		vsets.push_back(v_set);
		ng++;
	}
	cout<<"sets: "<<vsets.size()<<endl;
	//classifying
	vector<bool> tag(vsets.size(), true);
	vector<vector<int>> esets(vsets.size());
	for(auto &e: g_.Es)
	{
		if(v_flag[e.vs[0]] == v_flag[e.vs[1]])
		{
			int group_id = v_flag[e.vs[0]];
			esets[group_id].push_back(e.id);

			if(e.color != Edge_tag::R)
				tag[group_id] = false;
		}
	}
	//splitting
	vector<int> v_map(g_.Vs.size(), 0), rv_map(g_.Vs.size(), 0);
	int n = vsets.size();
	for(int i = 0;i< n;i++)
	{
		//if(tag[i])//if pure red
			continue;

		Graph gt;
		gt.dimension = g_.dimension;
		
		gt.Vs.resize(vsets[i].size());
		Graph_Node gn;		
		for(int j=0;j<vsets[i].size();j++)
		{
			gn.id = j;
			gn.boundary = g_.Vs[vsets[i][j]].boundary;
			gn.position = g_.Vs[vsets[i][j]].position;
			gn.normal = g_.Vs[vsets[i][j]].normal;
			gn.direction = g_.Vs[vsets[i][j]].direction;
			gn.scale = g_.Vs[vsets[i][j]].scale;
			gn.v = g_.Vs[vsets[i][j]].v;
			v_map[vsets[i][j]] = j;
			rv_map[j] = vsets[i][j];
			gt.Vs[j] = gn;
		}
		gt.Es.resize(esets[i].size());		
		Graph_Edge ge;
		ge.vs.resize(2);
		for(int j=0;j<esets[i].size();j++)
		{
			auto &e = g_.Es[esets[i][j]];
			ge.vs[0] = v_map[e.vs[0]];
			ge.vs[1] = v_map[e.vs[1]];
			ge.boundary = e.boundary;
			ge.id = j;
			ge.time = 0;
			ge.color = e.color;
			ge.w = e.w;
			ge.posy_dif = e.posy_dif;
			gt.Es[j] = ge;

			gt.Vs[ge.vs[0]].nvs.insert(ge.vs[1]);
			gt.Vs[ge.vs[1]].nvs.insert(ge.vs[0]);
			gt.Vs[ge.vs[0]].nes.insert(ge.id);
			gt.Vs[ge.vs[1]].nes.insert(ge.id);
		}
		vector<vector<int>> vsets_;
		//split_group(gt, vsets_);
		//split_group_based_on_color(gt, vsets_);

		if(vsets_.size()>1)
		{
			for(auto &vs: vsets_)for(auto &vid: vs) vid = rv_map[vid];
			vsets[i] = vsets_[0];
			vsets.insert(vsets.end(), vsets_.begin()+1, vsets_.end());
			for (int j = 1; j < vsets_.size(); j++)
				tag.push_back(false);
			tag[i] = false;

			cout <<"#s "<< vsets_.size() << endl;
		}
	}
	cout<<"sets: "<<vsets.size()<<endl;

	Graph gnew;
	gnew.dimension = g_.dimension;
	std::fill(v_map.begin(), v_map.end(), -1);
	for(int i=0;i<vsets.size();i++)
	{
		auto &vs = vsets[i];
		if(vs.size()<=1)
			continue;

		Graph_Node gn;
		gn.id = gnew.Vs.size();
		gn.singularity = !tag[i];
		group_representative(g_, vs, gn);
		gnew.Vs.push_back(gn);

		for(const auto &vid:vs)
			v_map[vid] = gn.id;
	}
	//extraction
	vector<tuple<int, int, int>> es;
	tuple<int, int, int> te;
	for(auto &e:g_.Es)
	{		
		int v0 = v_map[e.vs[0]];
		int v1 = v_map[e.vs[1]];
		if(v0 == v1)
			continue;
		if(v0 == -1 || v1 == -1)
			continue;
		if(v0 > v1)
			std::swap(v0, v1);

		es.push_back(make_tuple(e.id, v0, v1));
	}

	auto is_equal = [&](const std::tuple<int, int, int>& e0, const std::tuple<int, int, int>& e1) -> bool
	{
		return std::get<1>(e0) == std::get<1>(e1) && std::get<2>(e0) == std::get<2>(e1);
	};

	std::sort(es.begin(), es.end(), [](const std::tuple<int, int, int>& e0, const std::tuple<int, int, int>& e1) {
		return (std::get<1>(e0) < std::get<1>(e1)) || (std::get<1>(e0) == std::get<1>(e1) && std::get<2>(e0) < std::get<2>(e1));
	});


	Graph_Edge ge;
	ge.vs.resize(2);
	for(int i = 0; i < es.size(); i++)
	{
		auto &e = es[i];
		if(i==0 || !is_equal(es[i-1], es[i]))
		{
			int eid = get<0>(e);
			ge.id = gnew.Es.size();
			ge.boundary = false;
			ge.boundary = g_.Es[eid].boundary;
			ge.vs[0] = get<1>(e);
			ge.vs[1] = get<2>(e);
			gnew.Es.push_back(ge);

			gnew.Vs[ge.vs[0]].nvs.insert(ge.vs[1]);
			gnew.Vs[ge.vs[1]].nvs.insert(ge.vs[0]);
			gnew.Vs[ge.vs[0]].nes.insert(ge.id);
			gnew.Vs[ge.vs[1]].nes.insert(ge.id);
		}else 
			gnew.Es[ge.id].boundary = g_.Es[get<0>(e)].boundary || gnew.Es[ge.id].boundary;
	}
	g_ = gnew;
	cout<<"vs: "<<g_.Vs.size()<<" es "<<g_.Es.size()<<endl;
}
void MultiResolutionHierarchy::re_coloring_diagonals(Graph &g_)
{
	for (auto &v : g_.Vs)
	{
		//if (v.boundary) continue;

		vector<bool> tags;
		if (g_.dimension == 2)
			tags.resize(4);
		else if (g_.dimension == 3)
			tags.resize(6);

		MatrixXf D(tags.size(), 3);
		if (g_.dimension == 2)
		{
			Vector3f d = v.direction;
			D.row(0) = d;
			D.row(1) = d.cross(v.normal);
			D.row(2) = -D.row(0);
			D.row(3) = -D.row(1);
		}else if (g_.dimension == 3)
		{
			Quaternion q = v.direction;
			Eigen::Matrix<Float, 3, 3> M;
			M = q.toMatrix();
			D.block(0, 0, 3, 3) = M;
			D.block(3, 0, 3, 3) = -M;
		}
		//tagging blue edges
		std::fill(tags.begin(), tags.end(), false);
		for (auto &eid : v.nes)
		{
			auto &e = g_.Es[eid];
			if (e.color == Edge_tag::B)
			{
				Vector3f ed = g_.Vs[e.vs[0]].position - g_.Vs[e.vs[1]].position;
				if (e.vs[1] != v.id) ed *= -1;
				ed.normalize();

				VectorXf result = D *ed;
				VectorXf::Index id;
				result.maxCoeff(&id);
				tags[id] = true;
			}
		}
		//scan for not-captured direction
		bool all_captured = true;
		for (int i = 0; i < tags.size(); i++)
		{
			if (!tags[i]) {
				all_captured = false;
				break;
			}
		}
		if (all_captured)
			continue;
		//tagging non-blue edges
		vector<int> ntags(tags.size(),-1);
		vector<double> coefs(tags.size(), -1);
		for (auto &eid : v.nes)
		{
			auto &e = g_.Es[eid];
			if (e.color != Edge_tag::B)
			{
				Vector3f ed = g_.Vs[e.vs[0]].position - g_.Vs[e.vs[1]].position;
				if (e.vs[1] != v.id) ed *= -1;
				ed.normalize();

				VectorXf result = D *ed;
				VectorXf::Index id;
				result.maxCoeff(&id);
				double max_dot = result.maxCoeff();

				if ((ntags[id] == -1) || (ntags[id] != -1 && coefs[id] < max_dot))
				{
					ntags[id] = e.id;
					coefs[id] = max_dot;
				}
			}
		}

		if (g_.dimension == 3 && v.boundary)
		{
			int ben = 0;
			for (auto &eid : v.nes)
			{
				auto &e = g_.Es[eid];
				if (e.boundary && e.color == Edge_tag::B)
					ben++;
			}
			if (ben == 4) continue;
		}
		for (int i = 0; i < tags.size(); i++)
		{
			if (!tags[i] && ntags[i] != -1) {
				g_.Es[ntags[i]].color = Edge_tag::B;
				//cout << "eid: "<<ntags[i]<<" to blue: "<< g_.Es[ntags[i]].vs[0]<<" "<< g_.Es[ntags[i]].vs[1] << endl;
			}
		}
	}
}
void MultiResolutionHierarchy::remove_diagonals(Graph &g_)
{
	Graph gnew;
	gnew.dimension = g_.dimension;
	gnew.Es = g_.Es;
	g_.Es.clear();

	for (auto &v:g_.Vs)
	{
		v.nvs.clear();
		v.nes.clear();
	}
	Graph_Edge ge;
	for (auto &e: gnew.Es)
	{
		if ((g_.dimension == 2 && e.boundary) || e.color == Edge_tag::B)
		{
			ge = e;
			ge.id = g_.Es.size();
			g_.Es.push_back(ge);

			g_.Vs[ge.vs[0]].nvs.insert(ge.vs[1]);
			g_.Vs[ge.vs[1]].nvs.insert(ge.vs[0]);
			g_.Vs[ge.vs[0]].nes.insert(ge.id);
			g_.Vs[ge.vs[1]].nes.insert(ge.id);
		}
	}
	cout << "after removing diagonals ---- vs: " << g_.Vs.size() << " es " << g_.Es.size() << endl;

	gnew.Vs = g_.Vs;
	gnew.Es = g_.Es;
	g_.Vs.clear();
	g_.Es.clear();

	vector<int> v_map(gnew.Vs.size(), -1);	
	Graph_Node gn;
	for (auto &v : gnew.Vs)
	{
		if (v.nvs.size() == 1)// || (v.nvs.size() ==2 && gnew.dimension == 3))
			continue;

		gn = v;		
		gn.id = g_.Vs.size();		
		g_.Vs.push_back(gn);

		v_map[v.id] = gn.id;

		v.nvs.clear();
		v.nes.clear();
	}
	for (auto &e : gnew.Es)
	{
		int v0 = v_map[e.vs[0]];
		int v1 = v_map[e.vs[1]];
		if (v0 == -1 || v1 == -1)
			continue;
		if (v0 > v1)
			std::swap(v0, v1);

		ge = e;
		ge.vs[0] = v0;
		ge.vs[1] = v1;
		ge.id = g_.Es.size();
		g_.Es.push_back(ge);

		g_.Vs[ge.vs[0]].nvs.insert(ge.vs[1]);
		g_.Vs[ge.vs[1]].nvs.insert(ge.vs[0]);
		g_.Vs[ge.vs[0]].nes.insert(ge.id);
		g_.Vs[ge.vs[1]].nes.insert(ge.id);
	}
	cout << "after removing dangling edges ---- vs: " << g_.Vs.size() << " es " << g_.Es.size() << endl;

}
void MultiResolutionHierarchy::re_projection(Graph &g_)
{
	if (!D3) return;
	MatrixXu F(3, m_sur.Fs.size());
	for (auto &f : m_sur.Fs) {
		F(0, f.id) = f.vs[0];
		F(1, f.id) = f.vs[1];
		F(2, f.id) = f.vs[2];
	}
	BVH *bvh_local = new BVH(&F, &m_sur.V, mAABB);
	bvh_local->build();

	std::vector<int> V_flag(g_.Vs.size(), 0);
	for (auto &gn : g_.Vs)
	{
		if (!gn.boundary) continue;
		V_flag[gn.id] = 1;
		double scale = 3 * gn.scale[0];
		if (bvh_local && bvh_local->F()->size() > 0) {
			V_flag[gn.id] = 2;
			Ray ray1(gn.position, gn.normal, 0, scale);
			Ray ray2(gn.position, -gn.normal, 0, scale);
			uint32_t idx1 = 0, idx2 = 0;
			Float t1 = 0, t2 = 0;
			bvh_local->rayIntersect(ray1, idx1, t1);
			bvh_local->rayIntersect(ray2, idx2, t2);
			if (std::min(t1, t2) < scale)
			{
				V_flag[gn.id] = 3;
				gn.position = t1 < t2 ? ray1(t1) : ray2(t2);
			}
			else
			{
				;// cout << gn.position << " gn.normal " << gn.normal << " scale" << scale << endl;
			}
			
		}
	}
	//write_Vertex_Types_TXT(V_flag, "C:/xgao/meshing/tensor-field-graphing/data/3d/v_flag.txt");

}
void MultiResolutionHierarchy::graph_extraction_fast()
{
	bool changes = true;
	int iter = 0;
	edge_tagging(g);
	composit_edges_colors(g, E_rend);
	grouping_collapse(g);
	edge_tagging(g);
	re_coloring_diagonals(g);
	remove_diagonals(g);
	re_projection(g);
	// int vn_pre = g.Vs.size(), en_pre = g.Es.size();
	// while(changes)
	// {
	// 	std::cout<<iter++<<" iteration..."<<std::endl;
	// 	changes = false;


	// 	// connect_close_points(g);
	// 	// vector<int> ledges;
	// 	// detect_parallel_redundant_edges(g, ledges);
	// 	// split_parallel_redundant_edges(g, ledges);
	// 	// cout <<"split "<< ledges.size()/2 << endl;
		
	// 	if (vn_pre == g.Vs.size() && en_pre == g.Es.size())
	// 		changes = false;
	// 	else
	// 	{
	// 		vn_pre = g.Vs.size();
	// 		en_pre = g.Es.size();
	// 	}
	// }
	//remove_dangling(g);

	mV_final.resize(3, g.Vs.size());
	mQ_final.resize(3, g.Vs.size());
	mN_final.resize(3, g.Vs.size());

	for (auto &v : g.Vs) {
		mV_final.col(v.id) = v.position;
		mQ_final.col(v.id) = v.direction;
		mN_final.col(v.id) = v.normal;
	}
}

