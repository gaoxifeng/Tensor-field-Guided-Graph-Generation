#pragma once
#include "common.h"

extern bool mesh_preprocessing(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, VectorXb &BC, Mesh &m);
extern bool data_preprocessing(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, VectorXb &BC, bool & D3);
extern void write_surface_mesh_OBJ(Mesh &m, string &path);
extern void orient_polygon_mesh(MatrixXf &HV, std::vector<std::vector<uint32_t>> &HF);
extern void extract_surface_mesh(Mesh &meshi, Mesh &mesho);
extern bool load_Jun_2D_3D_nonUniform(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, bool & D3);
extern bool load_Jun_2D_3D(const std::string &path, MatrixXu &E, MatrixXf &V, VectorXb &BV, MatrixXf &N, MatrixXf &C, MatrixXf &Q, MatrixXf &S, bool & D3);
extern bool load_Jun_2D_3D(const std::string &path, MatrixXu &E, MatrixXf &V, MatrixXf &Q, MatrixXf &S, bool & D3);
extern bool write_Graph(const std::string &path,const MatrixXf &V, const MatrixXu &E);
extern bool write_Graph(const std::string &path, const Graph &g_);

extern bool load_Jun(const std::string &path, MatrixXu &F, MatrixXf &V, MatrixXf &Q, MatrixXf &S);
extern bool load_Jun(const std::string &path, MatrixXu &F, MatrixXf &V, MatrixXf &Q);

extern void load_obj(const std::string &filename, MatrixXu &F, MatrixXf &V);
extern void write_surface_mesh_OBJ(MatrixXf &V, MatrixXu &F, char * path);

extern void write_Vertex_Types_TXT(std::vector<int> &V_types, string path);

extern void write_statistics_TXT(statistics &sta, char * path);
