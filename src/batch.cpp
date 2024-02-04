/*
    batch.cpp -- command line interface to Instant Meshes

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#include "batch.h"
#include "meshio.h"
#include "hierarchy.h"
#include "optimizer_sequencial.h"
#include "timer.h"

void batch_process(char *input, char *output,
	uint32_t dimension, Float tlen, Float scale, int smooth_iter, bool rosy) {
	//cout << endl;
	//cout << "Running in batch mode:" << endl;
	//cout << "   Input file             = " << input << endl;
	//cout << "   Output file            = " << output << endl;
	char patho[1024];

	MultiResolutionHierarchy mRes;

	Timer<> timer;
	timer.beginStage("data pre-processing");
	mRes.load(input);
	
	mRes.ratio_scale = mRes.ms.mAverageEdgeLength * 3.0 * 1.0/scale / mRes.diagonalLen;

	timer.endStage();
	mRes.build();
     
	//mRes.setScale(scale);
	Optimizer_sequencial *mOptimizer;
	mOptimizer = new Optimizer_sequencial(mRes);

	if (rosy)
	{
		timer.beginStage("rosy optimization");

		mOptimizer->setExtrinsic(true);
		mOptimizer->setMaxIterations(smooth_iter);
		mOptimizer->setOptimizeOrientations(true);
		mOptimizer->run();
		timer.endStage();
	}

	timer.beginStage("posy optimization");

	mOptimizer->setOptimizePositions(true);
	mOptimizer->run();
	timer.endStage();
	
	timer.beginStage("mesh extraction");
	
	mRes.construct_Graph(mRes.m, mRes.g);
	mRes.graph_extraction_fast();

	timer.endStage();

	//sprintf(patho, "%s%s", input, "_graph.obj");
	//write_Graph(patho, mRes.mV_final, mRes.mE_final);
	sprintf(patho, "%s%s", input, "_graph_opt.obj");
	write_Graph(patho, mRes.g);

	//system("PAUSE");
}
