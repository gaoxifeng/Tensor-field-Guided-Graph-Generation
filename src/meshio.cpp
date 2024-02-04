#include "meshio.h"
#include <fstream>
#include <unordered_map>
#include "timer.h"
#include <iomanip>
#include "orientations.h"
using namespace std;

void write_surface_mesh_OBJ(Mesh &m, string &path)
{
	std::fstream f(path, std::ios::out);
	for (auto &v : m.Vs)
		f << "v " << v.v[0] << " " << v.v[1] << " " << v.v[2] << std::endl;
	if (m.Fs.size())
	{
		for (auto &ff : m.Fs)
		{
			f << "f";
			for (auto &vid : ff.vs)
				f << " " << vid + 1;
			f << endl;
		}
	}
	else if (m.Es.size())
	{
		for (auto &e : m.Es) {
			f << "l";
			for (auto &vid : e.vs)
				f << " " << vid + 1;
			f << endl;
		}
	}
	f.close();
}

void build_connectivity(Mesh &hmi)
{
	if (hmi.type == Mesh_type::Qua || hmi.type == Mesh_type::Tri) {
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
		temp.reserve(hmi.Fs.size() * 3);
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i) {
			int vn = hmi.Fs[i].vs.size();
			for (uint32_t j = 0; j < vn; ++j) {
				uint32_t v0 = hmi.Fs[i].vs[j], v1 = hmi.Fs[i].vs[(j + 1) % vn];
				if (v0 > v1) std::swap(v0, v1);
				temp.push_back(std::make_tuple(v0, v1, i, j));
			}
			hmi.Fs[i].es.resize(vn);
		}
		std::sort(temp.begin(), temp.end());
		hmi.Es.reserve(temp.size() / 2);
		uint32_t E_num = 0;
		Hybrid_E e; e.boundary = true; e.vs.resize(2);
		for (uint32_t i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				hmi.Es.push_back(e);
			}
			else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
				std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
				hmi.Es[E_num - 1].boundary = false;

			hmi.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : hmi.Vs) v.boundary = false;
		for (uint32_t i = 0; i < hmi.Es.size(); ++i)
			if (hmi.Es[i].boundary) {
				hmi.Vs[hmi.Es[i].vs[0]].boundary = hmi.Vs[hmi.Es[i].vs[1]].boundary = true;
			}
	}
	else if(hmi.type == Mesh_type::Tet)
	{
		std::vector<std::vector<uint32_t>> total_fs(hmi.Hs.size() * 4);
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> tempF(hmi.Hs.size() * 4);
		std::vector<uint32_t> vs(3);
		for (uint32_t i = 0; i < hmi.Hs.size(); ++i) {
			for (short j = 0; j < 4; j++) {
				for (short k = 0; k < 3; k++) vs[k] = hmi.Hs[i].vs[tet_faces[j][k]];
				uint32_t id = 4 * i + j;
				total_fs[id] = vs;
				std::sort(vs.begin(), vs.end());
				tempF[id] = std::make_tuple(vs[0], vs[1], vs[2], id, i, j);
			}
			hmi.Hs[i].fs.resize(4);
		}
		std::sort(tempF.begin(), tempF.end());
		hmi.Fs.reserve(tempF.size() / 3);
		Hybrid_F f; f.boundary = true;
		uint32_t F_num = 0;
		for (uint32_t i = 0; i < tempF.size(); ++i) {
			if (i == 0 || (i != 0 &&
				(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) || std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
					std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1])))) {
				f.id = F_num; F_num++;
				f.vs = total_fs[std::get<3>(tempF[i])];
				hmi.Fs.push_back(f);
			}
			else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) && std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
				std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1])))
				hmi.Fs[F_num - 1].boundary = false;

			hmi.Hs[std::get<4>(tempF[i])].fs[std::get<5>(tempF[i])] = F_num - 1;
		}

		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp(hmi.Fs.size() * 3);
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i) {
			for (uint32_t j = 0; j < 3; ++j) {
				uint32_t v0 = hmi.Fs[i].vs[j], v1 = hmi.Fs[i].vs[(j + 1) % 3];
				if (v0 > v1) std::swap(v0, v1);
				temp[3 * i + j] = std::make_tuple(v0, v1, i, j);
			}
			hmi.Fs[i].es.resize(3);
		}
		std::sort(temp.begin(), temp.end());
		hmi.Es.reserve(temp.size() / 2);
		uint32_t E_num = 0;
		Hybrid_E e; e.boundary = false; e.vs.resize(2);
		for (uint32_t i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				hmi.Es.push_back(e);
			}
			hmi.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : hmi.Vs) v.boundary = false;
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i)
			if (hmi.Fs[i].boundary) for (uint32_t j = 0; j < 3; ++j) {
				uint32_t eid = hmi.Fs[i].es[j];
				hmi.Es[eid].boundary = true;
				hmi.Vs[hmi.Es[eid].vs[0]].boundary = hmi.Vs[hmi.Es[eid].vs[1]].boundary = true;
			}
	}
	else if(hmi.type == Mesh_type::Hex) {

		std::vector<std::vector<uint32_t>> total_fs(hmi.Hs.size() * 6);
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> tempF(hmi.Hs.size() * 6);
		std::vector<uint32_t> vs(4);
		for (uint32_t i = 0; i < hmi.Hs.size(); ++i) {
			for (short j = 0; j < 6; j++){
				for (short k = 0; k < 4; k++) vs[k] = hmi.Hs[i].vs[hex_face_table[j][k]];
				uint32_t id = 6 * i + j;
				total_fs[id] = vs;
				std::sort(vs.begin(), vs.end());
				tempF[id] = std::make_tuple(vs[0], vs[1], vs[2], vs[3], id, i, j);
			}
			hmi.Hs[i].fs.resize(6);
		}
		std::sort(tempF.begin(), tempF.end());
		hmi.Fs.reserve(tempF.size() / 3);
		Hybrid_F f; f.boundary = true;
		uint32_t F_num = 0;
		for (uint32_t i = 0; i < tempF.size(); ++i) {
			if (i == 0 || (i != 0 &&
				(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) || std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
					std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1]) || std::get<3>(tempF[i]) != std::get<3>(tempF[i - 1])))) {
				f.id = F_num; F_num++;
				f.vs = total_fs[std::get<4>(tempF[i])];
				hmi.Fs.push_back(f);
			}
			else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) && std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
				std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1]) && std::get<3>(tempF[i]) == std::get<3>(tempF[i - 1])))
				hmi.Fs[F_num - 1].boundary = false;

			hmi.Hs[std::get<5>(tempF[i])].fs[std::get<6>(tempF[i])] = F_num - 1;
		}

		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp(hmi.Fs.size() * 4);
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i) {
			for (uint32_t j = 0; j < 4; ++j) {
				uint32_t v0 = hmi.Fs[i].vs[j], v1 = hmi.Fs[i].vs[(j + 1) % 4];
				if (v0 > v1) std::swap(v0, v1);
				temp[4 * i + j] = std::make_tuple(v0, v1, i, j);
			}
			hmi.Fs[i].es.resize(4);
		}
		std::sort(temp.begin(), temp.end());
		hmi.Es.reserve(temp.size() / 2);
		uint32_t E_num = 0;
		Hybrid_E e; e.boundary = false; e.vs.resize(2);
		for (uint32_t i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				hmi.Es.push_back(e);
			}
			hmi.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : hmi.Vs) v.boundary = false;
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i)
			if (hmi.Fs[i].boundary) for (uint32_t j = 0; j < 4; ++j) {
				uint32_t eid = hmi.Fs[i].es[j];
				hmi.Es[eid].boundary = true;
				hmi.Vs[hmi.Es[eid].vs[0]].boundary = hmi.Vs[hmi.Es[eid].vs[1]].boundary = true;
			}
	}
	//f_nhs;
	for (auto &f : hmi.Fs)f.neighbor_hs.clear();
	for (uint32_t i = 0; i < hmi.Hs.size(); i++) {
		for (uint32_t j = 0; j < hmi.Hs[i].fs.size(); j++) hmi.Fs[hmi.Hs[i].fs[j]].neighbor_hs.push_back(i);
	}
	//e_nfs, v_nfs
	for (auto &e : hmi.Es) e.neighbor_fs.clear();
	for (auto &v : hmi.Vs) v.neighbor_fs.clear();
	for (uint32_t i = 0; i < hmi.Fs.size(); i++) {
		for (uint32_t j = 0; j < hmi.Fs[i].vs.size(); j++) hmi.Vs[hmi.Fs[i].vs[j]].neighbor_fs.push_back(i);			
		for (uint32_t j = 0; j < hmi.Fs[i].es.size(); j++) hmi.Es[hmi.Fs[i].es[j]].neighbor_fs.push_back(i);
		 
	}
	//v_nes, v_nvs
	for (auto &v : hmi.Vs) {
		v.neighbor_es.clear();
		v.neighbor_vs.clear();
	}
	for (uint32_t i = 0; i < hmi.Es.size(); i++) {
		uint32_t v0 = hmi.Es[i].vs[0], v1 = hmi.Es[i].vs[1];
		hmi.Vs[v0].neighbor_es.push_back(i);
		hmi.Vs[v1].neighbor_es.push_back(i);
		hmi.Vs[v0].neighbor_vs.push_back(v1);
		hmi.Vs[v1].neighbor_vs.push_back(v0);
	}
	//e_nhs
	for (auto &e : hmi.Es) e.neighbor_hs.clear();
	for (uint32_t i = 0; i < hmi.Es.size(); i++) {
		std::vector<uint32_t> nhs;
		for (uint32_t j = 0; j < hmi.Es[i].neighbor_fs.size(); j++) {
			uint32_t nfid = hmi.Es[i].neighbor_fs[j];
			nhs.insert(nhs.end(), hmi.Fs[nfid].neighbor_hs.begin(), hmi.Fs[nfid].neighbor_hs.end());
		}
		std::sort(nhs.begin(), nhs.end()); nhs.erase(std::unique(nhs.begin(), nhs.end()), nhs.end());
		hmi.Es[i].neighbor_hs = nhs;
	}
	//v_nhs
	for (auto &v : hmi.Vs) v.neighbor_hs.clear();
	for (uint32_t i = 0; i < hmi.Hs.size(); i++) {
		for (uint32_t j = 0; j < hmi.Hs[i].vs.size(); j++) hmi.Vs[hmi.Hs[i].vs[j]].neighbor_hs.push_back(i);
	}
}

bool load_mesh(const std::string &path, Mesh &hmi)
{
	string path_mesh = path + ".obj";
	string path_info = path + ".txt";
	std::ifstream ff;
	ff.open(path_mesh, std::ios::in);
	if (ff.fail())
	{
		path_mesh = path + ".mesh";
		ff.open(path_mesh, std::ios::in);
		if (ff.fail()) 
			return false;
		else 
			hmi.type = Mesh_type::Tet;
	}
	else
		hmi.type = Mesh_type::Tri;

	std::ifstream f_info;
	f_info.open(path_info, std::ios::in);
	if(ff.fail())
		return false;

	std::cout<<"start to read"<<endl;
	hmi.Vs.clear();
	hmi.Es.clear();
	hmi.Fs.clear();
	hmi.Hs.clear();

	if (hmi.type == Mesh_type::Tri)
	{
		int vnum, hnum;	float x, y, z;
		char s[1024];
		ff.getline(s, 1023);
		std::sscanf(s, "%d %d", &vnum, &hnum);

		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			std::sscanf(s, "%f %f %f", &x, &y, &z);

			Hybrid_V v;
			v.id = i;
			v.boundary = false;
			v.v.push_back(x);
			v.v.push_back(y);
			v.v.push_back(z);			

			f_info.getline(s, 1023);
			sscanf(s, "%f %f", &x, &y);
			v.direction = Vector3f(x, y, 0);
			f_info.getline(s, 1023);
			f_info.getline(s, 1023);
			sscanf(s, "%f %f", &x, &y);
			v.scale = Vector3f(abs(x),abs(y), 0);
			
			v.hard_constraint = false;
			v.normal = Vector3f(0, 0, 1);

			hmi.Vs.push_back(v);	
		}

		Hybrid_F f;
		f.vs.resize(3);
		uint32_t a, b, c, d, e;
		hmi.Fs.resize(hnum);

		for(int i = 0; i < hnum; i++)
		{
			ff.getline(s, 1023);
			sscanf(s, "%d %d %d %d", &vnum, &a, &b, &c);

			f.vs[0]=a;
			f.vs[1]=b;
			f.vs[2]=c;
			f.id = i; hmi.Fs[i] = f;
			for (uint32_t i = 0; i < f.vs.size(); i++) hmi.Vs[f.vs[i]].neighbor_fs.push_back(f.id);
		}
	}
	else if (hmi.type == Mesh_type::Tet) 
	{
		int vnum, hnum;	float x, y, z;
		char s[1024];
		ff.getline(s, 1023);
		std::sscanf(s, "%d %d", &vnum, &hnum);

		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			std::sscanf(s, "%f %f %f", &x, &y, &z);

			Hybrid_V v;
			v.id = i;
			v.boundary = false;
			v.v.push_back(x);
			v.v.push_back(y);
			v.v.push_back(z);			

			Eigen::Matrix<Float, 3, 3> M;

			f_info.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			M.col(0) = Vector3f(x, y, z); 
			M.col(0).normalize();

			f_info.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			M.col(1) = Vector3f(x, y, z);
			M.col(1).normalize();

			f_info.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);

			M.col(2) = M.col(0).cross(M.col(1));
			v.direction = Quaternion::toQuaternion(M).normalized();

			f_info.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			v.scale = Vector3f(abs(x),abs(y), z);

			v.hard_constraint = false;
			v.normal = Vector3f(0, 0, 1);

			hmi.Vs.push_back(v);	
		}

		Hybrid h;
		h.vs.resize(4);
		uint32_t a, b, c, d;
		hmi.Hs.resize(hnum);

		for(int i = 0; i < hnum; i++)
		{
			ff.getline(s, 1023);
			sscanf(s, "%d %d %d %d %d", &vnum, &a, &b, &c, &d);

			h.vs[0]=a;
			h.vs[1]=b;
			h.vs[2]=c;
			h.vs[3]=d;
			h.id = i; hmi.Hs[i] = h;
			for (uint32_t i = 0; i < h.vs.size(); i++) hmi.Vs[h.vs[i]].neighbor_hs.push_back(h.id);
		}
	}
	else
		return false;

	ff.close();
	std::cout<<"end of the reading"<<endl;	
	return true;
}
void extract_surface_mesh(Mesh &meshi, Mesh &mesho) {
	mesho.Vs.clear(); mesho.Es.clear(); mesho.Fs.clear(); mesho.Hs.clear();
	mesho.type = Mesh_type::Tri;

	vector<bool> V_tag(meshi.Vs.size(), false);
	vector<int32_t> V_map(meshi.Vs.size(), -1), V_map_reverse;

	for (auto f : meshi.Fs) if (f.boundary) {
		for (auto vid : f.vs) V_tag[vid] = true;

		if (f.vs.size() == 3) {
			Hybrid_F hf; 
			hf.id = mesho.Fs.size();
			hf.vs = f.vs;
			mesho.Fs.push_back(hf);
		}
		else if (f.vs.size() == 4) {
			Hybrid_F hf;
			hf.id = mesho.Fs.size();

			hf.vs.push_back(f.vs[0]);
			hf.vs.push_back(f.vs[1]);
			hf.vs.push_back(f.vs[2]);
			mesho.Fs.push_back(hf);
			hf.vs.clear();
			hf.vs.push_back(f.vs[2]);
			hf.vs.push_back(f.vs[3]);
			hf.vs.push_back(f.vs[0]);
			mesho.Fs.push_back(hf);
		}

	}
	//re-indexing
	uint32_t newV_id = 0;
	for (uint32_t i = 0; i < V_tag.size(); i++) if (V_tag[i]) {
		V_map[i] = newV_id++; V_map_reverse.push_back(i);
	}
	mesho.V.resize(3, newV_id);
	for (uint32_t i = 0; i < V_tag.size(); i++) if (V_tag[i]) {
		Hybrid_V v;
		v.v.push_back(meshi.V(0, i));
		v.v.push_back(meshi.V(1, i));
		v.v.push_back(meshi.V(2, i));
		v.id = mesho.Vs.size(); mesho.Vs.push_back(v);
		mesho.V.col(v.id) = meshi.V.col(i);
	}
	for (uint32_t i = 0; i < mesho.Fs.size(); i++) for (uint32_t j = 0; j < 3; j++) mesho.Fs[i].vs[j] = V_map[mesho.Fs[i].vs[j]];
	//orient direction
	build_connectivity(mesho);
}

void orient_polygon_mesh(MatrixXf &HV, vector<vector<uint32_t>> &HF) {
	//orient surface
	if (!HF.size())return;

	vector<vector<uint32_t>> HFE(HF.size());
	vector<tuple_E> Es;

	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>> temp;
	temp.reserve(HF.size() * 3);
	HFE.resize(HF.size());
	for (uint32_t f = 0; f < HF.size(); ++f) {
		for (uint32_t e = 0; e < HF[f].size(); ++e) {
			uint32_t v0 = HF[f][e], v1 = HF[f][(e + 1) % HF[f].size()];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f, e, Edge_tag::B));
		}
		std::vector<uint32_t> fes(HF[f].size());
		HFE[f] = fes;
	}
	std::sort(temp.begin(), temp.end());
	Es.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			Es.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), true, 0, std::get<4>(temp[i]), E_num, -1, 0));
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			std::get<2>(Es[E_num]) = false;
		HFE[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
	vector<vector<uint32_t>> Efs(Es.size());
	for (uint32_t i = 0; i < HFE.size(); i++)for (auto e : HFE[i])Efs[e].push_back(i);

	uint32_t start_f = 0;
	vector<bool> flag(HF.size(), true);
	flag[start_f] = false;
	std::queue<uint32_t> pf_temp; pf_temp.push(start_f);
	while (!pf_temp.empty()) {
		uint32_t fid = pf_temp.front(); pf_temp.pop();
		for (auto eid : HFE[fid]) for (auto nfid : Efs[eid]) {
			if (!flag[nfid]) continue;
			uint32_t v0 = std::get<0>(Es[eid]), v1 = std::get<1>(Es[eid]);
			int32_t v0_pos = std::find(HF[fid].begin(), HF[fid].end(), v0) - HF[fid].begin();
			int32_t v1_pos = std::find(HF[fid].begin(), HF[fid].end(), v1) - HF[fid].begin();

			if ((v0_pos + 1) % HF[fid].size() != v1_pos) swap(v0, v1);

			int32_t v0_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v0) - HF[nfid].begin();
			int32_t v1_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v1) - HF[nfid].begin();

			if ((v0_pos_ + 1) % HF[nfid].size() == v1_pos_) std::reverse(HF[nfid].begin(), HF[nfid].end());

			pf_temp.push(nfid); flag[nfid] = false;
		}
		if (pf_temp.empty()) {
			bool found = false;
			for (uint32_t i = 0; i < flag.size(); i++)if (flag[i]) {
				start_f = i;
				flag[start_f] = false; pf_temp.push(start_f);
				found = true;
			}
			if (!found) break;
		}
	}

	Float res = 0;
	Vector3f ori; ori.setZero();
	for (uint32_t i = 0; i < HF.size(); i++) {
		auto &fvs = HF[i];
		Vector3f center; center.setZero(); for (auto vid : fvs) center += HV.col(vid); center /= fvs.size();

		for (uint32_t j = 0; j < fvs.size(); j++) {
			Vector3f x = HV.col(fvs[j]) - ori, y = HV.col(fvs[(j + 1) % fvs.size()]) - ori, z = center - ori;
			res += -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));
		}
	}
	if (res > 0) {
		for (uint32_t i = 0; i < HF.size(); i++) std::reverse(HF[i].begin(), HF[i].end());
	}
}
bool mesh_preprocessing(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, VectorXb &BC, Mesh &m)
{
	if(!load_mesh(path, m)) 
		return false;
	build_connectivity(m);

	m.V.resize(3, m.Vs.size());
	m.BV.resize(m.Vs.size());
	m.N.resize(3, m.Vs.size());
	m.C.resize(3, m.Vs.size());
	m.Q.resize(3, m.Vs.size());
	if(m.type == Mesh_type::Tet)
		m.Q.resize(4, m.Vs.size());
	m.S.resize(3, m.Vs.size());	
	m.BC.resize(m.Vs.size());
	
	for (const auto &v : m.Vs)
	{
		for (int j = 0; j < 3; j++)
			m.V(j, v.id) = v.v[j];
		m.BV[v.id] = v.boundary;
		m.N.col(v.id) = v.normal;
		if(m.type == Mesh_type::Tet)
			m.N.col(v.id) = Vector3f::Zero();
		m.Q.col(v.id) = v.direction;
		m.S.col(v.id) = v.scale;
		m.BC[v.id] = v.hard_constraint;
	}
	for (const auto &v : m.Vs)
	{
		if(m.type == Mesh_type::Tri){			
			m.C.col(v.id) = Eigen::Vector3f::Zero();
		
			if(!v.boundary) continue;

			std::vector<uint32_t> vs;
			for (auto eid : v.neighbor_es) {
				if (!m.Es[eid].boundary) continue;

				uint32_t v0 = m.Es[eid].vs[0];
				uint32_t v1 = m.Es[eid].vs[1];
				if (v0 == v.id) vs.push_back(v1);
				else vs.push_back(v0);
			}
			Vector3f direct0, direct1; direct1.setZero();
			direct0 = (m.V.col(v.id) - m.V.col(vs[0])).normalized();
			for (uint32_t k = 1; k < vs.size(); k++) {
				Vector3f direct_ = (m.V.col(vs[k]) - m.V.col(v.id)).normalized();			
				direct1 += direct_;
			}
			m.C.col(v.id) = (direct0 + direct1).normalized();	
			if(m.C.col(v.id) == Eigen::Vector3f::Zero())		
			m.C.col(v.id) = Eigen::Vector3f(1,0,0);
		}
	}
	if(m.type == Mesh_type::Tet)
	{
		int nbf = 0;
		for(auto &f:m.Fs)
			if(f.boundary) nbf++;

		vector<vector<uint32_t>> HF(nbf);
		MatrixXu mF(3, nbf); nbf = 0;
		for(auto &f:m.Fs)
			if(f.boundary){
				//mF.col(f_b) = Fs[f];
				HF[nbf].push_back(f.vs[0]);
				HF[nbf].push_back(f.vs[1]);
				HF[nbf].push_back(f.vs[2]);
				nbf++;
			}

		orient_polygon_mesh(m.V, HF);

		for (uint32_t i = 0; i < HF.size(); i++)
			for (uint32_t j = 0; j < 3; j++)mF(j, i) = HF[i][j];

		MatrixXf mNF, mCF;
		mNF.resize(3, mF.cols()); mCF.resize(3, mF.cols());
		VectorXi count(m.V.cols());
		count.setZero();
		for (uint32_t i = 0; i<mF.cols(); ++i) {
			uint32_t i0 = mF(0, i), i1 = mF(1, i), i2 = mF(2, i);
			Vector3f v0 = m.V.col(i0), v1 = m.V.col(i1), v2 = m.V.col(i2);
			Vector3f n = (v1 - v0).cross(v2 - v0).normalized();
			mNF.col(i) = n;
			mCF.col(i) += (v0 + v1 + v2) / 3;
			m.N.col(i0) += n; m.N.col(i1) += n; m.N.col(i2) += n;
			count[i0]++; count[i1]++; count[i2]++;
		}

		for (uint32_t i = 0; i<m.N.cols(); ++i) {
			m.Vs[i].normal = m.N.col(i).normalized();
			
			if(m.Vs[i].boundary&& (m.N.col(i) == Vector3f::Zero()|| std::isnan(m.N.col(i)[0])))
			{
				cout<<"normal error"<< count[i] <<" "<<m.N.col(i)<< endl;
				m.N.col(i) = m.Vs[i].normal = Vector3f(1, 0, 0);				
			}

			m.C.col(i) = Vector3f::Zero();
			if (m.N.col(i) != Vector3f::Zero()) {
				Vector3f d1 = m.N.col(i) / count[i],
					d2 = m.N.col(i).normalized();
				if (d1.dot(d2) < 0.85f)
					d2 = Vector3f::Zero();;//
				m.N.col(i) = d2;
				if(std::isnan(m.N.col(i)[0]))
				{
					
					cout<<"normal error"<<count[i]<<"; "<<m.Vs[i].normal<<"; "<<m.N.col(i)<< endl;
				}					
				if (d2 != Vector3f::Zero())
					m.C.col(i) = m.V.col(i);
			}
		}
	}
	BV = m.BV;
	C = m.C;
	Q = m.Q;
	N = m.N;
	V = m.V;
	S = m.S;	
	BC = m.BC;

	//S.setConstant(1);
	//for (auto &v : m.Vs)
	//	v.scale = S.col(v.id);
	//Float min = S.block(0, 0, 2, S.cols()).minCoeff(), max = S.block(0, 0, 2, S.cols()).maxCoeff();
	//for (int i = 0; i < S.rows(); i++)
	//	for (int j = 0; j < S.cols(); j++)
	//		S(i, j) = 1.0 / S(i, j) * max;
	//for (auto &v : m.Vs)
	//	v.scale = S.col(v.id);
	
	E.resize(2, m.Es.size());
	for(int i = 0; i<m.Es.size();i++){
		E(0, i) = m.Es[i].vs[0];
		E(1, i) = m.Es[i].vs[1];
	}

	return true;
}

void averaging_properties(Mesh &M_, std::vector<uint32_t> &vs, Hybrid_V &v)
{
	VectorXf q_i = M_.Vs[vs[0]].direction;
	Quaternion q_i_, q_j_;
	if(M_.type == Mesh_type::Hex)
	 	q_i_= q_i;

	Vector3f n_i = M_.Vs[vs[0]].normal;
	Vector3f s_i = M_.Vs[vs[0]].scale;
	float d_i = M_.Vs[vs[0]].density;
	for (int j = 1;j<vs.size();j++) {
		VectorXf q_j = M_.Vs[vs[j]].direction;
		if(M_.type == Mesh_type::Hex)
	 		q_j_= q_j;
		Vector3f n_j = M_.Vs[vs[j]].normal;
		Vector3f s_j = M_.Vs[vs[j]].scale;
		float d_j = M_.Vs[vs[j]].density;
		int which = 0;
		
		if(M_.type == Mesh_type::Hex)
		{
			auto q_j__ = Quaternion::applyRotation(q_i_, q_j_);
			Eigen::Matrix<Float, 3, 3> M0, M1;
			M0 = q_i_.toMatrix();
			M1 = q_j__.toMatrix();
			for(int k=0;k<3;k++){
				int id=0;Float dmax=-1;
				for(int m=0;m<3;m++){
					Float d = std::abs(M0.col(k).dot(M1.col(m)));
					if(d>dmax){dmax = d; id = m;}	
				}
				s_j[k] = M_.Vs[vs[j]].scale[id];
			}
			q_i += q_j__; 
			s_i += s_j;
		}else
		{
			which = findRotation(q_i, n_i, q_j, n_j);
			auto q_j_ = applyRotationKeep(q_i, n_i, q_j, n_j);
			q_i += q_j_ - n_i * n_i.dot(q_j_); 
			if(which==1||which ==3){
				s_i[0] += s_j[1];
				s_i[1] += s_j[0];
				s_i[2] += s_j[2];
			}else s_i += s_j;
		}
		d_i += d_j;

	}
	s_i /= vs.size();
	d_i /= vs.size();
	if ((M_.type == Mesh_type::Qua && q_i == Eigen::Vector3f::Zero())||
	(M_.type == Mesh_type::Hex && q_i_ == Quaternion::Zero())) {
		q_i = M_.Vs[vs[0]].direction;
		n_i = M_.Vs[vs[0]].normal;
		s_i = M_.Vs[vs[0]].scale;
	}else q_i.normalize();
	v.direction = q_i;
	v.scale = s_i;
	v.normal = n_i;
	v.density = d_i;
	v.hard_constraint = false;
}
void refine_input(Mesh &M, const int iter)
{
	if(iter <=0)
		return;

	for (int i = 0; i < iter; i++) {
		cout<<"iter "<<i<<endl;
		Mesh M_;
		M_.type = M.type;

		std::vector<int> E2V(M.Es.size()), F2V(M.Fs.size()), H2V(M.Hs.size());

		int vn = 0;
		for (auto v : M.Vs) {
			Hybrid_V v_;
			v_.id = vn++;
			v_.v.resize(3);
			for (int j = 0; j < 3; j++)v_.v[j] = M.V(j, v.id);

			v_.direction = v.direction;
			v_.scale = v.scale;
			v_.normal = v.normal;
			v_.density = v.density;
			v_.hard_constraint = v.hard_constraint;

			M_.Vs.push_back(v_);
		}
		
		for (auto e : M.Es) {
			Hybrid_V v;
			v.id = vn++;
			v.v.resize(3);

			Vector3f center;
			center.setZero();
			for (auto vid: e.vs) center += M.V.col(vid);
			center /= e.vs.size();
			for (int j = 0; j < 3; j++)v.v[j] = center[j];

			averaging_properties(M, e.vs, v);

			M_.Vs.push_back(v);
			E2V[e.id] = v.id;
		}
		for (auto f : M.Fs) {
			Hybrid_V v;
			v.id = vn++;
			v.v.resize(3);

			Vector3f center;
			center.setZero();
			for (auto vid : f.vs) center += M.V.col(vid);
			center /= f.vs.size();
			for (int j = 0; j < 3; j++)v.v[j] = center[j];
			
			if(f.vs.size() !=4)
			{
				std::cout<<"fid"<<f.id<<" "<< f.vs.size()<<endl;
				for(auto &vid: f.vs) cout<<vid<<endl;
			}
			averaging_properties(M, f.vs, v);

			M_.Vs.push_back(v);
			F2V[f.id] = v.id;
		}
		for (auto h : M.Hs) {
			Hybrid_V v;
			v.id = vn++;
			v.v.resize(3);

			Vector3f center;
			center.setZero();
			for (auto vid : h.vs) center += M.V.col(vid);
			center /= h.vs.size();
			for (int j = 0; j < 3; j++)v.v[j] = center[j];

			averaging_properties(M, h.vs, v);

			M_.Vs.push_back(v);
			H2V[h.id] = v.id;
		}
		if(M_.type == Mesh_type::Qua)
		{
			for(auto f : M.Fs) {
				auto v = M_.Vs[F2V[f.id]];

				Hybrid_F f_;
				f_.vs.resize(4);
				for(int j=0;j<4;j++)
				{
					f_.id = M_.Fs.size();
					f_.vs[0] = f.vs[j];

					auto A_ = M.Vs[f.vs[j]].neighbor_es;
					auto B_ = M.Vs[f.vs[(j+1)%4]].neighbor_es;
					std::vector<int> shared_es;
					std::sort(A_.begin(), A_.end());
					std::sort(B_.begin(), B_.end());
					std::set_intersection(A_.begin(), A_.end(), B_.begin(), B_.end(), std::back_inserter(shared_es));

					f_.vs[1] = E2V[shared_es[0]];
					f_.vs[2] = v.id;

					B_ = M.Vs[f.vs[(j+3)%4]].neighbor_es;
					std::sort(B_.begin(), B_.end());
					shared_es.clear();
					std::set_intersection(A_.begin(), A_.end(), B_.begin(), B_.end(), std::back_inserter(shared_es));

					f_.vs[3] = E2V[shared_es[0]];

					M_.Fs.push_back(f_);
				}
			}
		}
		else if(M_.type == Mesh_type::Hex)
		{
			for(auto h : M.Hs) {
				for (auto vid : h.vs) {
					//top 4 vs
					vector<uint32_t> top_vs(4);
					top_vs[0] = vid;
					int fid = -1;
					for (auto nfid : M.Vs[vid].neighbor_fs)if (find(h.fs.begin(), h.fs.end(), nfid) != h.fs.end()) {
						fid = nfid; break;
					}
					assert(fid != -1);
					top_vs[2] = F2V[fid];
					
					int v_ind = find(M.Fs[fid].vs.begin(), M.Fs[fid].vs.end(),vid) - M.Fs[fid].vs.begin();
					int e_pre = M.Fs[fid].es[(v_ind - 1 + 4) % 4];
					int e_aft = M.Fs[fid].es[v_ind];
					top_vs[1] = E2V[e_pre];
					top_vs[3] = E2V[e_aft];
					//bottom 4 vs
					vector<uint32_t> bottom_vs(4);

					auto nvs = M.Vs[vid].neighbor_vs;
					int e_per = -1;
					for (auto nvid : nvs) if (find(h.vs.begin(), h.vs.end(), nvid) != h.vs.end()) {
						if (nvid != M.Es[e_pre].vs[0] && nvid != M.Es[e_pre].vs[1] &&
							nvid != M.Es[e_aft].vs[0] && nvid != M.Es[e_aft].vs[1]) {
							vector<uint32_t> sharedes, es0 = M.Vs[vid].neighbor_es, es1 = M.Vs[nvid].neighbor_es;
							sort(es0.begin(), es0.end()); sort(es1.begin(), es1.end());
							set_intersection(es0.begin(), es0.end(), es1.begin(), es1.end(), back_inserter(sharedes));
							assert(sharedes.size());
							e_per = sharedes[0];
							break;
						}
					}

					assert(e_per != -1);
					bottom_vs[0] = E2V[e_per];
					bottom_vs[2] = H2V[h.id];

					int f_pre = -1;
					vector<uint32_t> sharedfs, fs0 = M.Es[e_pre].neighbor_fs, fs1 = M.Es[e_per].neighbor_fs;
					
					sort(fs0.begin(), fs0.end()); sort(fs1.begin(), fs1.end());
					set_intersection(fs0.begin(), fs0.end(), fs1.begin(), fs1.end(), back_inserter(sharedfs));
					for (auto sfid : sharedfs)if (find(h.fs.begin(), h.fs.end(), sfid) != h.fs.end()) {
						f_pre = sfid; break;
					}
					assert(f_pre != -1);

					int f_aft = -1;
					sharedfs.clear();
					fs0 = M.Es[e_aft].neighbor_fs;
					fs1 = M.Es[e_per].neighbor_fs;
					sort(fs0.begin(), fs0.end()); sort(fs1.begin(), fs1.end());
					set_intersection(fs0.begin(), fs0.end(), fs1.begin(), fs1.end(), back_inserter(sharedfs));
					for (auto sfid : sharedfs)if (find(h.fs.begin(), h.fs.end(), sfid) != h.fs.end()) {
						f_aft = sfid; break;
					}
					assert(f_aft != -1);

					bottom_vs[1] = F2V[f_pre];
					bottom_vs[3] = F2V[f_aft];

					vector<uint32_t> ele_vs = top_vs;
					ele_vs.insert(ele_vs.end(), bottom_vs.begin(), bottom_vs.end());
					
					//new ele
					Hybrid h_;
					h_.id = M_.Hs.size();
					h_.vs = ele_vs;
					M_.Hs.push_back(h_);
				}
			}
		}

		build_connectivity(M_);

		M_.V.resize(3, M_.Vs.size());
		M_.BV.resize(M_.Vs.size());
		M_.N.resize(3, M_.Vs.size());
		M_.C.resize(3, M_.Vs.size());
		M_.Q.resize(3, M_.Vs.size());
		if(M_.type == Mesh_type::Hex)
			M_.Q.resize(4, M_.Vs.size());
		M_.S.resize(3, M_.Vs.size());
		M_.D.resize(M_.Vs.size());
		M_.BC.resize(M_.Vs.size());
		for (auto v : M_.Vs)
		{
			for (int j = 0; j < 3; j++)M_.V(j, v.id) = v.v[j];
			M_.Q.col(v.id) = v.direction;
			M_.S.col(v.id) = v.scale;
			M_.BV[v.id] = v.boundary;
			M_.N.col(v.id) = v.normal;
			M_.D[v.id] = v.density;
			M_.BC[v.id] = v.hard_constraint;
		}
		M = M_;
	}
}

bool load_data(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, VectorXf &D,VectorXb &BC,
std::map<std::vector<int>, int> &Coord2inds,
 bool & D3, int &UU, int &VV, int &WW)
{
	string in_format = "";
	size_t last_slash_idx = path.rfind('.');
	in_format = path.substr(last_slash_idx);
	if (in_format != ".Jun") return false;

	std::ifstream ff(path, std::ios::in);
	if (ff.fail())return false;
	std::cout<<"start to read"<<endl;
	bool D2_ = false, D3_ = false;
	char s[1024];
	int vnum;	float x, y, z;
	ff.getline(s, 1023);

	D3=false;
	if (sscanf(s, "%d %d %d", &UU, &VV, &WW) == 3) D3_ = true;
	else if (sscanf(s, "%d %d", &UU, &VV) == 2) D2_ = true;
	if(!D2_ && !D3_) return false;
	D3=D3_;

	// vnum = UU*VV;
	// if(WW!=0) vnum*=WW;
	ff.getline(s, 1023);
	sscanf(s, "%d", &vnum);	
	V.resize(3, vnum);
	S.resize(3, vnum);
	BV.resize(vnum);
	BV.setZero();
	N.resize(3, vnum);
	C.resize(3, vnum);
	D.resize(vnum);
	BC.resize(vnum);

	if(D2_){
		Q.resize(3, vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			sscanf(s, "%f %f", &x, &y);
			V.col(i) = Vector3f(x, y, 0);
			std::vector<int> coord;
			coord.push_back(x);
			coord.push_back(y);
			coord.push_back(0);
			Coord2inds[coord] = i;
			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			Q.col(i) = Vector3f(x, y, 0);
			ff.getline(s, 1023);
			ff.getline(s, 1023);

			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			// ff.getline(s, 1023);
			// sscanf(s, "%f", &z);
			z = y;
			S.col(i) = Vector3f(abs(x),abs(y), abs(z));
			//S.col(i) = Vector3f(abs(y), abs(x),abs(z));
			//density
			ff.getline(s, 1023);
			sscanf(s, "%f", &z);
			D[i] = z;
			//hard constraint
			int h;
			ff.getline(s, 1023);
			sscanf(s, "%d", &h);
			BC[i] = h;			
		}
	}
	if(D3_){
		//Q.resize(3, 3 * vnum);
		Q.resize(4, vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			V.col(i) = Vector3f(x, y, z);	
			std::vector<int> coord;
			coord.push_back(x);
			coord.push_back(y);
			coord.push_back(z);
			Coord2inds[coord] = i;

			Eigen::Matrix<Float, 3, 3> M;

			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			// Q.col(3 * i) = Vector3f(x, y, z);
			// Q.col(3 * i).normalize();
			M.col(0) = Vector3f(x, y, z); 

			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			// Q.col(3 * i + 1) = Vector3f(x, y, z);
			// Q.col(3 * i + 1).normalize();
			M.col(1) = Vector3f(x, y, z);

			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			// Q.col(3 * i + 2) = Vector3f(x, y, z);
			// Q.col(3 * i + 2).normalize();

			M.col(2) = M.col(0).cross(M.col(1));
			Q.col(i) = Quaternion::toQuaternion(M).normalized();


			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			//S.col(i) = Vector3f(abs(x),abs(y), abs(z));
			S.col(i) = Vector3f(1, 1, 1);

			ff.getline(s, 1023);
			sscanf(s, "%f", &z);
			D[i] = z;
			//hard constraint
			int h;
			ff.getline(s, 1023);
			sscanf(s, "%d", &h);
			BC[i] = h;			
		}
	}

	ff.close();
	std::cout<<"end of the reading"<<endl;	
	return true;
}
void density_removal(Mesh &mesh, double density_threshold)
{
	std::vector<bool> v_flag(mesh.Vs.size(), true), f_flag(mesh.Fs.size(), true), h_flag(mesh.Hs.size(), true);	
	std::vector<int> v_map(mesh.Vs.size(), -1);
	
	Mesh m_;
	m_.type =mesh.type;
	int vnum_ = 0;
	for(int i=0;i< mesh.Vs.size();i++)
	{
		if(mesh.D[i] < density_threshold)
			v_flag[i] = false;	
		else
		{
			v_map[i] =vnum_;
			vnum_++;

			Hybrid_V v;
			v.id = m_.Vs.size();
			v.boundary = false;
			v.v.push_back(mesh.V(0, i));
			v.v.push_back(mesh.V(1, i));
			v.v.push_back(mesh.V(2, i));
			v.direction = mesh.Q.col(i);
			v.scale = mesh.S.col(i);
			v.normal = mesh.N.col(i);
			v.density = mesh.D[i];
			v.hard_constraint = mesh.BC[i];
			m_.Vs.push_back(v);

		}	
	}

	if(mesh.type == Mesh_type::Qua)
	{
		for(auto f: mesh.Fs)
		{
			for(auto &vid: f.vs)	
				if(!v_flag[vid]) f_flag[f.id] =false;

			if(!f_flag[f.id]) continue;

			Hybrid_F f_;
			f_.id = m_.Fs.size();
			for(auto vid: f.vs)
				f_.vs.push_back(v_map[vid]);			
			m_.Fs.push_back(f_);
		}
	}
	else if(mesh.type == Mesh_type::Hex)
	{
		for(auto h: mesh.Hs)
		{
			for(auto &vid: h.vs)	
				if(!v_flag[vid]) h_flag[h.id] =false;

			if(!h_flag[h.id]) continue;

			Hybrid h_;
			h_.id = m_.Hs.size();
			for(auto vid: h.vs)
				h_.vs.push_back(v_map[vid]);			
			m_.Hs.push_back(h_);
		}
	}
	build_connectivity(m_);

	m_.V.resize(3, m_.Vs.size());
	m_.BV.resize(m_.Vs.size());
	m_.N.resize(3, m_.Vs.size());
	m_.C.resize(3, m_.Vs.size());
	m_.Q.resize(3, m_.Vs.size());
	if(m_.type == Mesh_type::Hex)
		m_.Q.resize(4, m_.Vs.size());
	m_.S.resize(3, m_.Vs.size());
	m_.D.resize(m_.Vs.size());
	m_.BC.resize(m_.Vs.size());
	for (auto v : m_.Vs)
	{
		for (int j = 0; j < 3; j++)m_.V(j, v.id) = v.v[j];
		m_.Q.col(v.id) = v.direction;
		m_.S.col(v.id) = v.scale;
		m_.BV[v.id] = v.boundary;
		m_.N.col(v.id) = v.normal;
		m_.D[v.id] = v.density;
		m_.BC[v.id] = v.hard_constraint;
	}

	mesh = m_;
}
void smoothing(Mesh &mesh, int iter)
{
	MatrixXi E;
	MatrixXf V = mesh.V.transpose();
	VectorXi B;
	int ne =0;
	for(auto &e: mesh.Es) if(e.boundary) ne++;
	E.resize(ne, 2); ne=0;
	for(auto &e: mesh.Es) if(e.boundary) {
		E(ne, 0) = e.vs[0];
		E(ne, 1) = e.vs[1];
		ne++;
	}
	int nv =0;
	for(auto &v: mesh.Vs) if(v.hard_constraint) nv++;
	B.resize(nv); nv = 0;
	for(auto &v: mesh.Vs) if(v.hard_constraint) B[nv++] = v.id;

	//smooth_::boundarySmoothing_2D(V, E, B);
	mesh.V = V.transpose();
	// //smoothing
	// for(int k = 0; k < iter; k++)
	// {
	// 	if(mesh.type == Mesh_type::Qua){
	// 		for(int i = 0; i < mesh.Vs.size(); i++)
	// 		{		
	// 			if(!mesh.Vs[i].boundary || mesh.Vs[i].hard_constraint) continue;
	// 			Vector3f centroid = mesh.V.col(i);
	// 			std::vector<uint32_t> vs;
	// 			for (auto eid : mesh.Vs[i].neighbor_es) {
	// 				if (!mesh.Es[eid].boundary) continue;

	// 				uint32_t v0 = mesh.Es[eid].vs[0];
	// 				uint32_t v1 = mesh.Es[eid].vs[1];
	// 				if (v0 == i) vs.push_back(v1);
	// 				else vs.push_back(v0);
	// 			}
	// 			for(auto &vid: vs)
	// 				centroid += mesh.V.col(vid);
	// 			if(vs.size()) centroid /= (1+vs.size());
	// 			mesh.V.col(i) = centroid;
	// 		}
	// 	}
	// }
}
bool data_preprocessing(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, VectorXb &BC, bool & D3)
{
	Mesh mesh;
	int UU = 0, VV = 0, WW = 0;
	if(!load_data(path, E, mesh.V, mesh.BV, mesh.N, mesh.C, mesh.Q, mesh.S, mesh.D, mesh.BC, mesh.Coord2inds, D3, UU, VV, WW))
		return false;
	mesh.type = Mesh_type::Qua;
	if(D3) mesh.type = Mesh_type::Hex;

	//quad-grid
	for(int i = 0;i<mesh.V.cols();i++)
	{
		Hybrid_V v;
		v.id = i;
		v.boundary = false;
		v.hard_constraint = mesh.BC[i];
		v.v.push_back(mesh.V(0, i));
		v.v.push_back(mesh.V(1, i));
		v.v.push_back(mesh.V(2, i));
		v.direction = mesh.Q.col(i);
		v.scale = mesh.S.col(i);
		v.normal = Vector3f(0, 0, 1);
		v.density = mesh.D[i];
		mesh.Vs.push_back(v);
	}

	Hybrid_F f;	
	f.vs.resize(4);
	Hybrid h;	
	h.vs.resize(8);
	std::vector<int> coord(3);
	for(int i = 0; i < mesh.Vs.size(); i++)
	{
		int a = (int)mesh.V(0, i) - 1, b = (int)mesh.V(1, i) - 1, c = (int)mesh.V(2, i) - 1;
		if(D3){
			h.id = mesh.Hs.size();
			a++, b++, c++;
			int a_, b_, c_, id;
			h.vs[0] = i;
			bool not_found= false;
			for(int j = 1; j< 8; j++)
			{
				switch(j)
				{
					case 1: a_ = a + 1, b_ = b    , c_ = c    ; break;
					case 2: a_ = a + 1, b_ = b + 1, c_ = c    ; break;
					case 3: a_ = a    , b_ = b + 1, c_ = c    ; break;
					case 4: a_ = a    , b_ = b    , c_ = c + 1; break;
					case 5: a_ = a + 1, b_ = b    , c_ = c + 1; break;
					case 6: a_ = a + 1, b_ = b + 1, c_ = c + 1; break;
					case 7: a_ = a    , b_ = b + 1, c_ = c + 1; break;
				}
				if(a_ <= 0 || b_ <= 0 || c_ <= 0) {not_found = true; break;}
				if(a_ > UU || b_ > VV || c_ > WW) {not_found = true; break;}
			
				coord[0] = a_;
				coord[1] = b_;
				coord[2] = c_;
				id = mesh.Coord2inds.find(coord)->second;
				h.vs[j] = id;
			}
			if(!not_found)
				mesh.Hs.push_back(h);
		}
		else
		{
			f.id = mesh.Fs.size();
			f.vs[0] = i;

			int a_, b_, id;
			a_ = a + 1, b_ = b;
			if(a_ < 0 || b_ < 0) continue;
			if(a_ >= UU || b_ >= VV) continue;			
			id = b_ * UU + a_;			
			f.vs[1] = id;

			a_ = a + 1, b_ = b + 1;
			if(a_ < 0 || b_ < 0) continue;
			if(a_ >= UU || b_ >= VV) continue;
			id = b_ * UU + a_;
			f.vs[2] = id;
			
			a_ = a, b_ = b + 1;
			if(a_ < 0 || b_ < 0) continue;
			if(a_ >= UU || b_ >= VV) continue;
			id = b_ * UU + a_;
			f.vs[3] = id;	
			mesh.Fs.push_back(f);
		}
	}
	build_connectivity(mesh);
	mesh.V.resize(3, mesh.Vs.size());
	mesh.BV.resize(mesh.Vs.size());
	mesh.N.resize(3, mesh.Vs.size());
	mesh.C.resize(3, mesh.Vs.size());
	mesh.Q.resize(3, mesh.Vs.size());
	if(mesh.type == Mesh_type::Hex)
		mesh.Q.resize(4, mesh.Vs.size());
	mesh.S.resize(3, mesh.Vs.size());
	mesh.D.resize(mesh.Vs.size());
	mesh.BC.resize(mesh.Vs.size());
	for (auto v : mesh.Vs)
	{
		for (int j = 0; j < 3; j++)mesh.V(j, v.id) = v.v[j];
		mesh.Q.col(v.id) = v.direction;
		mesh.S.col(v.id) = v.scale;
		mesh.BV[v.id] = v.boundary;
		mesh.N.col(v.id) = v.normal;
		mesh.D[v.id] = v.density;
		mesh.BC[v.id] = v.hard_constraint;
	}

	refine_input(mesh,0);

	float density_threshold = 0.5;	
	density_removal(mesh, density_threshold);
	//smoothing(mesh, 2);

	//C
	for(int i = 0; i < mesh.Vs.size(); i++){
		if(mesh.type == Mesh_type::Qua){			
			mesh.C.col(i) = Eigen::Vector3f::Zero();
			
			if(!mesh.Vs[i].boundary) continue;

			std::vector<uint32_t> vs;
			for (auto eid : mesh.Vs[i].neighbor_es) {
				if (!mesh.Es[eid].boundary) continue;

				uint32_t v0 = mesh.Es[eid].vs[0];
				uint32_t v1 = mesh.Es[eid].vs[1];
				if (v0 == i) vs.push_back(v1);
				else vs.push_back(v0);
			}
			Vector3f direct0, direct1; direct1.setZero();
			direct0 = (mesh.V.col(i) - mesh.V.col(vs[0])).normalized();
			for (uint32_t k = 1; k < vs.size(); k++) {
				Vector3f direct_ = (mesh.V.col(vs[k]) - mesh.V.col(i)).normalized();
				direct1 += direct_;
			}
			mesh.C.col(i) = (direct0 + direct1).normalized();			
		}
	}

	BV = mesh.BV;
	C = mesh.C;
	Q = mesh.Q;
	N = mesh.N;
	V = mesh.V;
	S = mesh.S;
	S.setConstant(1);
	BC = mesh.BC;

	E.resize(2, mesh.Es.size());
	for(int i = 0; i<mesh.Es.size();i++){
		E(0, i) = mesh.Es[i].vs[0];
		E(1, i) = mesh.Es[i].vs[1];
	}

	return true;
}
bool load_Jun_2D_3D_nonUniform(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, bool & D3) {
	string in_format = "";
	size_t last_slash_idx = path.rfind('.');
	in_format = path.substr(last_slash_idx);
	if (in_format != ".Jun") return false;

	std::ifstream ff(path, std::ios::in);
	if (ff.fail())return false;
	std::cout<<"start to read"<<endl;
	bool D2_ = false, D3_ = false;
	char s[1024];
	int UU=0, VV=0, WW = 0; int vnum;	float x, y, z;
	ff.getline(s, 1023);

	D3=false;
	if (sscanf(s, "%d %d %d", &UU, &VV, &WW) == 3) D3_ = true;
	else if (sscanf(s, "%d %d", &UU, &VV) == 2) D2_ = true;
	if(!D2_ && !D3_) return false;
	D3=D3_;

	vnum = UU*VV;
	if(WW!=0) vnum*=WW;
	V.resize(3, vnum);
	S.resize(3, vnum);
	BV.resize(vnum);
	BV.setZero();
	N.resize(3, vnum);
	C.resize(3, vnum);
	VectorXf D(vnum);

	Mesh mesh;
	if(D2_){
		Q.resize(3, vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			sscanf(s, "%f %f", &x, &y);
			V.col(i) = Vector3f(x, y, 0);	
			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			Q.col(i) = Vector3f(x, y, 0);
			ff.getline(s, 1023);
			ff.getline(s, 1023);

			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			// ff.getline(s, 1023);
			// sscanf(s, "%f", &z);
			z = y;
			S.col(i) = Vector3f(abs(x),abs(y), abs(z));
			//density
			ff.getline(s, 1023);
			sscanf(s, "%f", &z);
			D[i] = z;

			Hybrid_V v;
			v.id = i;
			v.boundary = false;
			mesh.Vs.push_back(v);
		}
	}
	if(D3_){
		ff.getline(s, 1023);
		Q.resize(3, 3 * vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			V.col(i) = Vector3f(x, y, z);	
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i) = Vector3f(x, y, z);
			Q.col(3 * i).normalize();
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i + 1) = Vector3f(x, y, z);
			Q.col(3 * i + 1).normalize();
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i + 2) = Vector3f(x, y, z);
			Q.col(3 * i + 2).normalize();

			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			S.col(i) = Vector3f(abs(x),abs(y), abs(z));
		}
	}

	float density_threshold = 0.5;
	
	//quad-grid
	Hybrid_F f;
	f.vs.resize(4);
	for(int i = 0; i < vnum; i++)
	{
		if(D[i]<density_threshold) continue;

		int a = (int)V(0, i) - 1, b = (int)V(1, i) - 1, c = (int)V(2, i);
		if(WW>0) c--;
		if(D3){
			a++;b++;c++;
		}

		if(D2_)
		{
			f.id = mesh.Fs.size();
			f.vs[0] = i;

			int a_, b_, id;
			a_ = a + 1, b_ = b;
			if(a_ < 0 || b_ < 0) continue;
			if(a_ >= UU || b_ >= VV) continue;			
			id = b_ * UU + a_;
			if(D[id]<density_threshold) continue;
			f.vs[1] = id;

			a_ = a + 1, b_ = b + 1;
			if(a_ < 0 || b_ < 0) continue;
			if(a_ >= UU || b_ >= VV) continue;
			id = b_ * UU + a_;
			if(D[id]<density_threshold) continue;
			f.vs[2] = id;
			
			a_ = a, b_ = b + 1;
			if(a_ < 0 || b_ < 0) continue;
			if(a_ >= UU || b_ >= VV) continue;
			id = b_ * UU + a_;
			if(D[id]<density_threshold) continue;
			f.vs[3] = id;	
		}

		for(auto &vid: f.vs)
			mesh.Vs[vid].neighbor_fs.push_back(f.id);
		mesh.Fs.push_back(f);
	}
	//extract edges
	if(D2_)
	{
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
		temp.reserve(mesh.Fs.size() * 4);
		for (uint32_t i = 0; i < mesh.Fs.size(); ++i) {
			int vn = mesh.Fs[i].vs.size();
			for (uint32_t j = 0; j < vn; ++j) {
				uint32_t v0 = mesh.Fs[i].vs[j], v1 = mesh.Fs[i].vs[(j + 1) % vn];
				if (v0 > v1) std::swap(v0, v1);
				temp.push_back(std::make_tuple(v0, v1, i, j));
			}
			mesh.Fs[i].es.resize(vn);
		}
		std::sort(temp.begin(), temp.end());
		mesh.Es.reserve(temp.size() / 2);
		uint32_t E_num = 0;
		Hybrid_E e; e.boundary = true; e.vs.resize(2);
		for (uint32_t i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				mesh.Es.push_back(e);
			}
			else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
				std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
				mesh.Es[E_num - 1].boundary = false;

			mesh.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : mesh.Vs) v.boundary = false;
		for (uint32_t i = 0; i < mesh.Es.size(); ++i)
		{
			if (mesh.Es[i].boundary) {
				mesh.Vs[mesh.Es[i].vs[0]].boundary = mesh.Vs[mesh.Es[i].vs[1]].boundary = true;
			}
			mesh.Vs[mesh.Es[i].vs[0]].neighbor_es.push_back(i);
			mesh.Vs[mesh.Es[i].vs[1]].neighbor_es.push_back(i);
		}
	}
	//smoothing
	for(int k = 0; k < 2; k++)
	{
		if(D2_){
			for(int i = 0; i < vnum; i++)
			{		
				if(!mesh.Vs[i].boundary) continue;
				Vector3f centroid = V.col(i);
				std::vector<uint32_t> vs;
				for (auto eid : mesh.Vs[i].neighbor_es) {
					if (!mesh.Es[eid].boundary) continue;

					uint32_t v0 = mesh.Es[eid].vs[0];
					uint32_t v1 = mesh.Es[eid].vs[1];
					if (v0 == i) vs.push_back(v1);
					else vs.push_back(v0);
				}
				for(auto &vid: vs)
					centroid += V.col(vid);
				if(vs.size()) centroid /= (1+vs.size());
				V.col(i) = centroid;
			}
		}
	}
	//BV, C, N
	for(int i = 0; i < vnum; i++){
		if(D2_){
			N.col(i) = Vector3f(0, 0, 1);
			BV[i] = mesh.Vs[i].boundary;
			C.col(i) = Eigen::Vector3f::Zero();
			
			if(!mesh.Vs[i].boundary) continue;

			std::vector<uint32_t> vs;
			for (auto eid : mesh.Vs[i].neighbor_es) {
				if (!mesh.Es[eid].boundary) continue;

				uint32_t v0 = mesh.Es[eid].vs[0];
				uint32_t v1 = mesh.Es[eid].vs[1];
				if (v0 == i) vs.push_back(v1);
				else vs.push_back(v0);
			}
			Vector3f direct0, direct1; direct1.setZero();
			direct0 = (V.col(i) - V.col(vs[0])).normalized();
			for (uint32_t k = 1; k < vs.size(); k++) {
				Vector3f direct_ = (V.col(vs[k]) - V.col(i)).normalized();
				direct1 += direct_;
			}
			C.col(i) = (direct0 + direct1).normalized();			
		}
	}
	//edges
	std::vector<std::vector<int>> Edges;
	std::vector<int> e(2);
	for(auto &v: mesh.Vs)
	{
		for(auto &fid: v.neighbor_fs)
		{
			for(auto &nvid: mesh.Fs[fid].vs)
			{
				if(nvid== v.id) continue;
				e[0] = v.id;
				e[1] = nvid;
				if(e[0]>e[1])std::swap(e[0], e[1]);

				Edges.push_back(e);
			}
		}
	}

	std::sort(Edges.begin(),Edges.end());
	Edges.erase(std::unique(Edges.begin(), Edges.end()), Edges.end());
	
	
	std::vector<bool> v_flag(vnum, true);
	std::vector<int> v_map(vnum, -1);
	int vnum_ = 0;
	for(int i=0;i<vnum;i++)
		if(D[i] < density_threshold)
			v_flag[i] = false;			
		else
		{
			v_map[i] =vnum_;
			vnum_++;
		}	
			
	auto V_ = V;
	auto S_ =S;
	auto BV_ =BV;
	auto N_ =N;
	auto C_ = C;
	auto Q_ = Q;
	V_.resize(3, vnum_);
	S_.resize(3, vnum_);
	BV_.resize(vnum_);
	BV_.setZero();
	N_.resize(3, vnum_);
	C_.resize(3, vnum_);
	Q_.resize(3, vnum_);
	for(int i=0;i<V.cols();i++)
	{
		if(v_map[i] != -1)
		{
			V_.col(v_map[i]) = V.col(i);
			S_.col(v_map[i]) = S.col(i);
			BV_[v_map[i]] = BV[i];
			N_.col(v_map[i]) = N.col(i);
			C_.col(v_map[i]) = C.col(i);
			Q_.col(v_map[i]) = Q.col(i);
		}
	}
	std::vector<std::vector<int>> Edges_;
	for(auto e: Edges)
	{
		if(v_map[e[0]]!=-1 && v_map[e[1]]!=-1)
		{
			e[0] = v_map[e[0]];
			e[1] = v_map[e[1]];
			Edges_.push_back(e);
		}
	}

	V = V_;
	S =S_;
	BV =BV_;
	N =N_;
	C = C_;
	Q = Q_;
	Edges = Edges_;

	E.resize(2, Edges.size());
	for(int i = 0; i<Edges.size();i++){
		E(0, i) = Edges[i][0];
		E(1, i) = Edges[i][1];
	}

	ff.close();
	std::cout<<"end of the reading"<<endl;	
	return true;
}


bool load_Jun_2D_3D(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, bool & D3) {


	string in_format = "";
	size_t last_slash_idx = path.rfind('.');
	in_format = path.substr(last_slash_idx);
	if (in_format != ".Jun") return false;

	std::ifstream ff(path, std::ios::in);
	if (ff.fail())return false;
	std::cout<<"start to read"<<endl;
	bool D2_ = false, D3_ = false;
	char s[1024];
	int UU=0, VV=0, WW = 0; int vnum;	float x, y, z;
	ff.getline(s, 1023);

	D3=false;
	if (sscanf(s, "%d %d %d", &UU, &VV, &WW) == 3) D3_ = true;
	else if (sscanf(s, "%d %d", &UU, &VV) == 2) D2_ = true;
	if(!D2_ && !D3_) return false;
	D3=D3_;

	vnum = UU*VV;
	if(WW!=0) vnum*=WW;
	V.resize(3, vnum);
	S.resize(3, vnum);
	BV.resize(vnum);
	BV.setZero();
	N.resize(3, vnum);
	C.resize(3, vnum);
	VectorXf D(vnum);

	if(D2_){
		Q.resize(3, vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			sscanf(s, "%f %f", &x, &y);
			V.col(i) = Vector3f(x, y, 0);	
			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			Q.col(i) = Vector3f(x, y, 0);
			ff.getline(s, 1023);
			ff.getline(s, 1023);

			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			ff.getline(s, 1023);
			sscanf(s, "%f", &z);
			S.col(i) = Vector3f(abs(x),abs(y), abs(z));
			//density
			ff.getline(s, 1023);
			sscanf(s, "%f", &z);
			D[i] = z;
		}
	}
	if(D3_){
		ff.getline(s, 1023);
		Q.resize(3, 3 * vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			V.col(i) = Vector3f(x, y, z);	
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i) = Vector3f(x, y, z);
			Q.col(3 * i).normalize();
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i + 1) = Vector3f(x, y, z);
			Q.col(3 * i + 1).normalize();
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i + 2) = Vector3f(x, y, z);
			Q.col(3 * i + 2).normalize();

			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			S.col(i) = Vector3f(abs(x),abs(y), abs(z));
		}
	}
	//edges
	std::vector<std::vector<int>> Edges;
	std::vector<int> e(2);
	for(int i = 0; i < vnum; i++){
		int a = (int)V(0, i) - 1, b = (int)V(1, i) - 1, c = (int)V(2, i);
		if(WW>0) c--;
		if(D3){
			a++;b++;c++;
		}

	//BV, N, C,
		if(D2_){
			N.col(i) = Vector3f(0, 0, 1);
			if((a==0 && b==0) || (a + 1 == UU && b + 1 == VV))
				C.col(i) = 	Vector3f(-std::sqrt(2)*0.5,std::sqrt(2)*0.5, 0);
			else if((a==0 && b + 1 == VV) || (a + 1 == UU && b==0))
				C.col(i) = 	Vector3f(std::sqrt(2)*0.5,std::sqrt(2)*0.5, 0);
			else if(a == 0 || a+ 1 == UU)
				C.col(i) = 	Vector3f(0,1, 0);
			else if(b == 0 || b+ 1 == VV)
				C.col(i) = 	Vector3f(1,0, 0);
			else 
				C.col(i) = 	Eigen::Vector3f::Zero();

			if(a == 0 || b == 0 || a + 1 == UU || b + 1 == VV)
				BV[i] = true;
		}
		else if(D3_){
			Vector3f n, dn;
			n.setZero(); 
			if(a==0)
				n += Vector3f(-1, 0, 0);
			else if(a +1 == UU)
				n += Vector3f(1, 0, 0);
			if(b==0)
				n += Vector3f(0, -1, 0);
			else if(b +1 == VV)
				n += Vector3f(0, 1, 0);
			if(c==0)
				n += Vector3f(0, 0, -1);
			else if(c +1 == WW)
				n += Vector3f(0, 0, 1);
			
			C.col(i) = Vector3f::Zero();
			if(n != Vector3f::Zero())
			{
				dn = n;
				n.normalize();
				if(abs(dn.dot(n))>0.85)
					C.col(i) = V.col(i);
			}
			N.col(i) = n;

			if(a == 0 || b == 0 || a + 1 == UU || b + 1 == VV || c ==0 || c+1 == WW)
				BV[i] = true;
		}
		

		int a_, b_, c_;
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				for(int l=0;l<3;l++){
					if(j==k && k==l && j==1) continue;
					if(WW==0 && l!=1) continue;
					a_ = a + j - 1;
					b_ = b + k - 1;
					c_ = c + l - 1;
					if(a_ < 0 || b_ < 0 || c_ < 0) continue;
					if(a_ >= UU || b_ >= VV || (c_ >= WW && WW!=0)) continue;
					
					int id = c_ * UU * VV + b_ * UU + a_;
					e[1] = id;
					e[0] = i;
					if(e[0]>e[1])std::swap(e[0],e[1]);

					Edges.push_back(e);
				}

	}
	std::sort(Edges.begin(),Edges.end());
	Edges.erase(std::unique(Edges.begin(), Edges.end()), Edges.end());
	
	float density_threshold = 0.5;
	std::vector<bool> v_flag(vnum, true);
	std::vector<int> v_map(vnum, -1);
	int vnum_ = 0;
	for(int i=0;i<vnum;i++)
		if(D[i] < density_threshold)
			v_flag[i] = false;			
		else
		{
			v_map[i] =vnum_;
			vnum_++;
		}	
			
	auto V_ = V;
	auto S_ =S;
	auto BV_ =BV;
	auto N_ =N;
	auto C_ = C;
	auto Q_ = Q;
	V_.resize(3, vnum_);
	S_.resize(3, vnum_);
	BV_.resize(vnum_);
	BV_.setZero();
	N_.resize(3, vnum_);
	C_.resize(3, vnum_);
	Q_.resize(3, vnum_);
	for(int i=0;i<V.cols();i++)
	{
		if(v_map[i] != -1)
		{
			V_.col(v_map[i]) = V.col(i);
			S_.col(v_map[i]) = S.col(i);
			BV_[v_map[i]] = BV[i];
			N_.col(v_map[i]) = N.col(i);
			C_.col(v_map[i]) = C.col(i);
			Q_.col(v_map[i]) = Q.col(i);
		}
	}
	std::vector<std::vector<int>> Edges_;
	for(auto e: Edges)
	{
		if(v_map[e[0]]!=-1 && v_map[e[1]]!=-1)
		{
			e[0] = v_map[e[0]];
			e[1] = v_map[e[1]];
			Edges_.push_back(e);
		}
	}

	V = V_;
	S =S_;
	BV =BV_;
	N =N_;
	C = C_;
	Q = Q_;
	Edges = Edges_;

	E.resize(2, Edges.size());
	for(int i = 0; i<Edges.size();i++){
		E(0, i) = Edges[i][0];
		E(1, i) = Edges[i][1];
	}

	ff.close();
	std::cout<<"end of the reading"<<endl;	
	return true;
}

bool load_Jun_2D_3D(const std::string &path, MatrixXu &E, MatrixXf &V, MatrixXf &Q, MatrixXf &S, bool & D3) {

	string in_format = "";
	size_t last_slash_idx = path.rfind('.');
	in_format = path.substr(last_slash_idx);
	if (in_format != ".Jun") return false;

	std::ifstream ff(path, std::ios::in);
	if (ff.fail())return false;
	std::cout<<"start to read"<<endl;
	bool D2_ = false, D3_ = false;
	char s[1024];
	int UU=0, VV=0, WW = 0; int vnum;	float x, y, z;
	ff.getline(s, 1023);

	D3=false;
	if (sscanf(s, "%d %d %d", &UU, &VV, &WW) == 3) D3_ = true;
	else if (sscanf(s, "%d %d", &UU, &VV) == 2) D2_ = true;
	if(!D2_ && !D3_) return false;
	D3=D3_;

	vnum = UU*VV;
	if(WW!=0) vnum*=WW;
	V.resize(3, vnum);
	S.resize(3, vnum);

	if(D2_){
		Q.resize(3, vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			sscanf(s, "%f %f", &x, &y);
			V.col(i) = Vector3f(x, y, 0);	
			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			Q.col(i) = Vector3f(x, y, 0);
			ff.getline(s, 1023);
			ff.getline(s, 1023);

			ff.getline(s, 1023);
			sscanf(s, "%f", &x);
			ff.getline(s, 1023);
			sscanf(s, "%f", &y);
			ff.getline(s, 1023);
			sscanf(s, "%f", &z);
			S.col(i) = Vector3f(abs(x),abs(y), abs(z));
		}
	}
	if(D3_){
		ff.getline(s, 1023);
		Q.resize(3, 3 * vnum);
		for (int i = 0; i < vnum; i++) {
			ff.getline(s, 1023);
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			V.col(i) = Vector3f(x, y, z);	
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i) = Vector3f(x, y, z);
			Q.col(3 * i).normalize();
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i + 1) = Vector3f(x, y, z);
			Q.col(3 * i + 1).normalize();
			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			Q.col(3 * i + 2) = Vector3f(x, y, z);
			Q.col(3 * i + 2).normalize();

			ff.getline(s, 1023);
			sscanf(s, "%f %f %f", &x, &y, &z);
			S.col(i) = Vector3f(abs(x),abs(y), abs(z));
		}
	}
	//edges
	std::vector<std::vector<int>> Edges;
	std::vector<int> e(2);
	for(int i = 0; i < vnum; i++){
		int a = (int)V(0, i) - 1, b = (int)V(1, i) - 1, c = (int)V(2, i);
		if(WW>0) c--;
		if(D3){
			a++;b++;c++;
		}
		
		int a_, b_, c_;
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				for(int l=0;l<3;l++){
					if(j==k && k==l && j==1) continue;
					if(WW==0 && l!=1) continue;
					a_ = a + j - 1;
					b_ = b + k - 1;
					c_ = c + l - 1;
					if(a_ < 0 || b_ < 0 || c_ < 0) continue;
					if(a_ >= UU || b_ >= VV || (c_ >= WW && WW!=0)) continue;
					
					int id = c_ * UU * VV + b_ * UU + a_;
					e[1] = id;
					e[0] = i;
					if(e[0]>e[1])std::swap(e[0],e[1]);

					Edges.push_back(e);
				}

	}
	std::sort(Edges.begin(),Edges.end());
	Edges.erase(std::unique(Edges.begin(), Edges.end()), Edges.end());
	
	E.resize(2, Edges.size());
	for(int i = 0; i<Edges.size();i++){
		E(0, i) = Edges[i][0];
		E(1, i) = Edges[i][1];
	}


	ff.close();
	std::cout<<"end of the reading"<<endl;	
	return true;
}
bool write_Graph(const std::string &path,const MatrixXf &V, const MatrixXu &E){
	std::fstream f(path, std::ios::out);
	//f << V.cols() << " "<<E.cols()<<endl; 
	for (int i = 0; i<V.cols(); i++)
		f << "v " << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	for (int i = 0; i < E.cols(); i++) {
		f << "l";
		for (int j = 0; j < E.rows(); j++)
			f << " " << E(j, i) + 1;
		f << endl;
	}
	f.close();
	return true;
}

bool write_Graph(const std::string &path, const Graph &g_) {
	std::fstream f(path, std::ios::out);
	//f << g_.Vs.size() << " " << g_.Es.size() << endl;
	for (int i = 0; i<g_.Vs.size(); i++)
		f << "v " << g_.Vs[i].position[0] << " " << g_.Vs[i].position[1] << " " << g_.Vs[i].position[2] << std::endl;

	for (int i = 0; i < g_.Es.size(); i++) {
		f << "l";
		for (int j = 0; j < g_.Es[i].vs.size(); j++)
			f << " " << g_.Es[i].vs[j] + 1;
		f << endl;
	}
	f.close();
	return true;
}

bool load_Jun(const std::string &path, MatrixXu &F, MatrixXf &V, MatrixXf &Q, MatrixXf &S) {

	string in_format = "";
	size_t last_slash_idx = path.rfind('.');
	in_format = path.substr(last_slash_idx);
	if (in_format != ".Jun") return false;

	std::ifstream ff(path, std::ios::in);
	if (ff.fail())return false;
	std::cout<<"start to read"<<endl;
	char s[1024];
	int r, c; int vnum, fnum;	float x, y, z;
	ff.getline(s, 1023);
	sscanf(s, "%d %d", &r, &c);

	vnum = r*c;
	fnum = 2 * (r-1)*(c-1);
	V.resize(3, vnum);
	F.resize(3, fnum);
	Q.resize(3, vnum);
	S.resize(3, vnum);

	for (int i = 0; i < vnum; i++) {
		ff.getline(s, 1023);
		sscanf(s, "%f %f", &x, &y);
		V.col(i) = Vector3f(x, y, 0);	
		ff.getline(s, 1023);
		sscanf(s, "%f", &x);
		ff.getline(s, 1023);
		sscanf(s, "%f", &y);
		Q.col(i) = Vector3f(x, y, 0);
		ff.getline(s, 1023);
		ff.getline(s, 1023);

		ff.getline(s, 1023);
		sscanf(s, "%f", &x);
		ff.getline(s, 1023);
		sscanf(s, "%f", &y);
		ff.getline(s, 1023);
		sscanf(s, "%f", &z);
		S.col(i) = Vector3f(abs(x),abs(y), abs(z));
	}
	
	fnum=0;
	for (int i = 1; i < c; i++)
		for (int j = 1; j < r; j++) {
			int v0= (i -1) * r + j-1, v1 = (i -1) * r + j, v2 = i * r + j;
			F.col(fnum++) = Vector3u(v0, v1, v2);
			v0= (i -1) * r + j-1, v1 = i * r + j, v2 = i * r + j - 1;
			F.col(fnum++) = Vector3u(v0, v1, v2);
		}

	ff.close();
	std::cout<<"end of the reading"<<endl;	
	return true;
}
bool load_Jun(const std::string &path, MatrixXu &F, MatrixXf &V, MatrixXf &Q) {

	string in_format = "";
	size_t last_slash_idx = path.rfind('.');
	in_format = path.substr(last_slash_idx);
	if (in_format != ".Jun") return false;

	std::ifstream ff(path, std::ios::in);
	if (ff.fail())return false;
	std::cout<<"start to read"<<endl;
	char s[1024];
	int r, c; int vnum, fnum;	float x, y, z;
	ff.getline(s, 1023);
	sscanf(s, "%d %d", &r, &c);

	vnum = r*c;
	fnum = 2 * (r-1)*(c-1);
	V.resize(3, vnum);
	F.resize(3, fnum);
	Q.resize(3, vnum);
	for (int i = 0; i < vnum; i++) {
		ff.getline(s, 1023);
		sscanf(s, "%f %f", &x, &y);
		V.col(i) = Vector3f(x, y, 0);	
		ff.getline(s, 1023);
		sscanf(s, "%f", &x);
		ff.getline(s, 1023);
		sscanf(s, "%f", &y);
		Q.col(i) = Vector3f(x, y, 0);
		ff.getline(s, 1023);
		ff.getline(s, 1023);
	}
	
	fnum=0;
	for (int i = 1; i < c; i++)
		for (int j = 1; j < r; j++) {
			int v0= (i -1) * r + j-1, v1 = (i -1) * r + j, v2 = i * r + j;
			F.col(fnum++) = Vector3u(v0, v1, v2);
			v0= (i -1) * r + j-1, v1 = i * r + j, v2 = i * r + j - 1;
			F.col(fnum++) = Vector3u(v0, v1, v2);
		}

	ff.close();
	std::cout<<"end of the reading"<<endl;	
	return true;
}
void load_obj(const std::string &filename, MatrixXu &F, MatrixXf &V) {
	/// Vertex indices used by the OBJ format
	struct obj_vertex {
		uint32_t p = (uint32_t)-1;
		uint32_t n = (uint32_t)-1;
		uint32_t uv = (uint32_t)-1;

		inline obj_vertex() { }

		inline obj_vertex(const std::string &string) {
			std::vector<std::string> tokens = str_tokenize(string, '/', true);

			if (tokens.size() < 1 || tokens.size() > 3)
				throw std::runtime_error("Invalid vertex data: \"" + string + "\"");

			p = str_to_uint32_t(tokens[0]);

#if 0
			if (tokens.size() >= 2 && !tokens[1].empty())
				uv = str_to_uint32_t(tokens[1]);

			if (tokens.size() >= 3 && !tokens[2].empty())
				n = str_to_uint32_t(tokens[2]);
#endif
		}

		inline bool operator==(const obj_vertex &v) const {
			return v.p == p && v.n == n && v.uv == uv;
		}
	};

	/// Hash function for obj_vertex
	struct obj_vertexHash : std::unary_function<obj_vertex, size_t> {
		std::size_t operator()(const obj_vertex &v) const {
			size_t hash = std::hash<uint32_t>()(v.p);
			hash = hash * 37 + std::hash<uint32_t>()(v.uv);
			hash = hash * 37 + std::hash<uint32_t>()(v.n);
			return hash;
		}
	};

	typedef std::unordered_map<obj_vertex, uint32_t, obj_vertexHash> VertexMap;

	size_t last_slash_idx = filename.rfind('.');
	if (filename.substr(last_slash_idx) != ".OBJ" && filename.substr(last_slash_idx) != ".obj") 
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");

	std::ifstream is(filename);
	if (is.fail())
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");
	cout << "Loading \"" << filename << "\" .. ";
	cout.flush();
	Timer<> timer;

	std::vector<Vector3f>   positions;
	//std::vector<Vector2f>   texcoords;
	//std::vector<Vector3f>   normals;
	std::vector<uint32_t>   indices;
	std::vector<obj_vertex> vertices;
	VertexMap vertexMap;

	std::string line_str;
	while (std::getline(is, line_str)) {
		std::istringstream line(line_str);

		std::string prefix;
		line >> prefix;

		if (prefix == "v") {
			Vector3f p;
			line >> p.x() >> p.y() >> p.z();
			positions.push_back(p);
		}
		else if (prefix == "vt") {
			/*
			Vector2f tc;
			line >> tc.x() >> tc.y();
			texcoords.push_back(tc);
			*/
		}
		else if (prefix == "vn") {
			/*
			Vector3f n;
			line >> n.x() >> n.y() >> n.z();
			normals.push_back(n);
			*/
		}
		else if (prefix == "f") {
			std::string v1, v2, v3, v4;
			line >> v1 >> v2 >> v3 >> v4;
			obj_vertex tri[6];
			int nVertices = 3;

			tri[0] = obj_vertex(v1);
			tri[1] = obj_vertex(v2);
			tri[2] = obj_vertex(v3);

			if (!v4.empty()) {
				/* This is a quad, split into two triangles */
				tri[3] = obj_vertex(v4);
				tri[4] = tri[0];
				tri[5] = tri[2];
				nVertices = 6;
			}
			/* Convert to an indexed vertex list */
			for (int i = 0; i<nVertices; ++i) {
				const obj_vertex &v = tri[i];
				VertexMap::const_iterator it = vertexMap.find(v);
				if (it == vertexMap.end()) {
					vertexMap[v] = (uint32_t)vertices.size();
					indices.push_back((uint32_t)vertices.size());
					vertices.push_back(v);
				}
				else {
					indices.push_back(it->second);
				}
			}
		}
	}
	F.resize(3, indices.size() / 3);
	memcpy(F.data(), indices.data(), sizeof(uint32_t)*indices.size());
	V.resize(3, vertices.size());
	for (uint32_t i = 0; i<vertices.size(); ++i)
		V.col(i) = positions.at(vertices[i].p - 1);
}
void write_surface_mesh_OBJ(MatrixXf &V, MatrixXu &F, char * path)
{
	std::fstream f(path, std::ios::out);
	for (int i = 0; i<V.cols(); i++)
		f << "v " << V(0, i) << " " << V(1, i) << " " << V(2, i) << std::endl;

	for (int i = 0; i < F.cols(); i++) {
		f << "f";
		for (int j = 0; j < F.rows(); j++)
			f << " " << F(j, i) + 1;
		f << endl;
	}
	f.close();
}
void write_Vertex_Types_TXT(std::vector<int> &V_types, string path) {
	std::fstream ff(path, std::ios::out);
	ff << V_types.size() << std::endl;
	for (auto type: V_types) {
			ff << std::fixed << type << endl;
	}
	ff.close();
}
void write_statistics_TXT(statistics &sta, char * path) {
	std::fstream f(path, std::ios::out);

	f<< 1 + 1 + 1 + sta.timings.size()+sta.polyhedral_ratios.size() << std::endl;
	f << sta.hex_ratio << "\t hex - ratio"<<std::endl;
	f << sta.tN << " " << sta.tetN << "\t #triangle #tet" << std::endl;
	f << sta.hN<<" "<<sta.pN <<"\t #hex #total"<< std::endl;
	//for(auto timing:sta.timings)
	if (sta.timings.size()) {
		f << sta.timings[0] / 1000 << "\t\t " << "timing: data pre-processing" << std::endl;
		f << sta.timings[1] / 1000 << "\t\t " << "timing: rosy optimization" << std::endl;
		f << sta.timings[2] / 1000 << "\t\t " << "timing: posy optimization" << std::endl;
		f << sta.timings[3] / 1000 << "\t\t " << "timing: mesh extraction" << std::endl;
	}

	for (uint32_t i = 0; i < sta.polyhedral_nums.size(); i++) {
		if(sta.polyhedral_nums[i]>0)
		f << sta.polyhedral_ratios[i] << "\t "<< i <<"th polyhedral "<<sta.polyhedral_nums[i] << std::endl;
	}
	f.close();
}
