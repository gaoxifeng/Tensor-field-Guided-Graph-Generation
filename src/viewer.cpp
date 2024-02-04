#include "viewer.h"
#include "im_resources.h"
#include "timer.h"
#include "bvh.h"
#include <nanogui/serializer/opengl.h>

Viewer::Viewer(std::string &filename, bool fullscreen)
    : Screen(Vector2i(1280, 960), "Robust Field-aligned Graph Generation", true),
      mOptimizer(nullptr) {

    mOptimizer = new Optimizer(mRes);

    /* Initialize shaders for rendering geometry and fields */
    mOrientationFieldShaderTet.init("orientation_field_shader_tet",
        (const char *) shader_orientation_field_tet_vert,
        (const char *) shader_orientation_field_tet_frag,
        (const char *) shader_orientation_field_tet_geo);

    mOrientationFieldShaderTri.init("orientation_field_shader_tri",
        (const char *) shader_orientation_field_tri_vert,
        (const char *) shader_orientation_field_tri_frag,
        (const char *) shader_orientation_field_tri_geo);

    mPositionFieldShader.init("position_field_shader",
        (const char *) shader_position_field_vert,
        (const char *) shader_position_field_frag);

    mGraphShader.init("mGraph_input",
        (const char *) shader_singularity_tet_vert,
        (const char *) shader_singularity_tet_frag);
	mEdge_color_morphShader.init("mEdge_color_morph",
		(const char *)shader_singularity_tet_vert,
			(const char *)shader_singularity_tet_frag);
	mEdge_color_morphRankedShader.init("mEdge_color_morphRanked",
		(const char *)shader_singularity_tet_vert,
		(const char *)shader_singularity_tet_frag);
	mExtractionResultShader.init("extraction_result",
        (const char *) shader_singularity_tet_vert,
        (const char *) shader_singularity_tet_frag);

    /* Default view setup */
    mCamera.arcball = Arcball();
    mCamera.arcball.setSize(mSize);
    mCamera.modelTranslation = -mRes.aabb().center().cast<float>();
    mCamera.modelZoom = 3.0f / mRes.aabb().extents().cwiseAbs().maxCoeff();
    mCamera.zoom = 1.0f;
    mLightPosition = Vector3f(0.0f, 0.3f, 5.0f);
    mBaseColor = Vector4f(0.4f, 0.5f, 0.7f, 1.f);
    mBaseColorBoundary = mRes.D3 ? Vector4f(0.0f, 0.0f, 1.0f, .2f) : mBaseColor;
    mSpecularColor = Vector4f(1.f, 1.f, 1.f, 1.f);
    mSpecularColorBoundary = mRes.D3 ? Vector4f(1.f, 1.f, 1.f, .2f) : mSpecularColor;
    mTranslate = false;

	/* Scan over example files in the 'datasets' directory */
	auto ctx = nvgContext();
	try {
		mExampleImages = nanogui::loadImageDirectory(ctx, "resources");
	}
	catch (const std::runtime_error &e) {
	}
	mExampleImages.insert(mExampleImages.begin(),
		std::make_pair(nvgImageIcon(ctx, loadmesh), ""));

	/* Initialize user interface */
    Window *window = new Window(this, "");
    window->setPosition(Vector2i(15, 15));
    window->setLayout(new GroupLayout());

	PopupButton *openBtn0 = new PopupButton(window, "Open mesh");
	openBtn0->setBackgroundColor(Color(0, 255, 0, 25));
	openBtn0->setIcon(ENTYPO_ICON_FOLDER);
	Popup *popup0 = openBtn0->popup();
	VScrollPanel *vscroll = new VScrollPanel(popup0);
	ImagePanel *panel0 = new ImagePanel(vscroll);
	panel0->setImages(mExampleImages);
	panel0->setCallback([&, openBtn0](int i) {
		openBtn0->setPushed(false);

		std::string filename2 = mExampleImages[i].second;
		std::string extension;
		if (filename2.size() > 4)
			extension = str_tolower(filename2.substr(filename2.size() - 4));

		if (filename2.empty()) {

			filename2 = nanogui::file_dialog({
				{ "obj", "Wavefront OBJ" },
				{ "vtk", "Visualization Toolkit VTk" },
				{ "mesh", "MESH" },
				{ "Jun", "general polyhedral mesh" }
			}, false);
			if (filename2 == "")
				return;
		}
		string path = "";
		size_t last_slash_idx = filename2.rfind('.');
		path = filename2.substr(0, last_slash_idx);

		if(!mRes.load(path)) return;
		mScaleBox->setValue(mRes.scale());

		filename = path;
		mRes.outpath = path;
		/* Default view setup */
		mCamera.arcball = Arcball();
		mCamera.arcball.setSize(mSize);
		mCamera.modelTranslation = -mRes.aabb().center().cast<float>();
		mCamera.modelZoom = 3.0f / mRes.aabb().extents().cwiseAbs().maxCoeff();
		mCamera.zoom = 1.0f;
		mLightPosition = Vector3f(0.0f, 0.3f, 5.0f);
		mBaseColor = Vector4f(0.4f, 0.5f, 0.7f, 1.f);
		//mBaseColor = Vector4f(1.0f, 1.0f, 1.0f, 1.f);
		mBaseColorBoundary = mRes.D3 ? Vector4f(0.0f, 0.0f, 1.0f, .2f) : mBaseColor;
		mSpecularColor = Vector4f(1.f, 1.f, 1.f, 1.f);
		mSpecularColorBoundary = mRes.D3 ? Vector4f(1.f, 1.f, 1.f, .2f) : mSpecularColor;
		mTranslate = false;

	});

	if (mRes.mV[0].cols()) {
		/* Default view setup */
		mCamera.arcball = Arcball();
		mCamera.arcball.setSize(mSize);
		mCamera.modelTranslation = -mRes.aabb().center().cast<float>();
		mCamera.modelZoom = 3.0f / mRes.aabb().extents().cwiseAbs().maxCoeff();
		mCamera.zoom = 1.0f;
		mLightPosition = Vector3f(0.0f, 0.3f, 5.0f);
		mBaseColor = Vector4f(0.4f, 0.5f, 0.7f, 1.f);
		//mBaseColor = Vector4f(1.0f, 1.0f, 1.0f, 1.f);
		mBaseColorBoundary = mRes.D3 ? Vector4f(0.0f, 0.0f, 1.0f, .2f) : mBaseColor;
		mSpecularColor = Vector4f(1.f, 1.f, 1.f, 1.f);
		mSpecularColorBoundary = mRes.D3 ? Vector4f(1.f, 1.f, 1.f, .2f) : mSpecularColor;
		mTranslate = false;
	}

//Parameters
	new Label(window, "Output scale", "sans-bold");
	mScaleBox = new FloatBox<Float>(window);
	mScaleBox->setValue(mRes.scale());
	mScaleBox->setEditable(true);
	mScaleBox->setAlignment(TextBox::Alignment::Right);
	mScaleBox->setId("outputscale");



//2D&3D
	Widget *statePanel = new Widget(window);
	statePanel->setLayout(new BoxLayout(Orientation::Horizontal, nanogui::Alignment::Middle, 0, 5));
	mSolveDatastructureBtn = new Button(statePanel, "Surface", ENTYPO_ICON_FLASH);
	mSolveDatastructureBtn->setBackgroundColor(Color(0, 0, 255, 25));
	mSolveDatastructureBtn->setFlags(Button::Flags::ToggleButton);
	mSolveDatastructureBtn->setChangeCallback([&](bool value) {
		mRes.build();

		auto const &R = mRes.E_input_rend;
		mGraphShader.bind();
		mGraphShader.uploadAttrib("position", MatrixXf(R.block(0, 0, 3, R.cols())));
		mGraphShader.uploadAttrib("color", MatrixXf(R.block(3, 0, 3, R.cols())));

		if (mRes.D3)
		{
			mOrientationFieldShaderTet.bind();
			mOrientationFieldShaderTet.uploadAttrib("position", mRes.V());
			mOrientationFieldShaderTet.uploadAttrib("q", mRes.Q());
			mOrientationFieldShaderTet.uploadAttrib("scale_value", mRes.mS[0]);
		}
		else
		{
			mOrientationFieldShaderTri.bind();
			mOrientationFieldShaderTri.uploadAttrib("position", mRes.V());
			mOrientationFieldShaderTri.uploadAttrib("q", mRes.Q());
			mOrientationFieldShaderTri.uploadAttrib("n", mRes.N());
			mOrientationFieldShaderTri.uploadAttrib("scale_value", mRes.mS[0]);
		}

		mPositionFieldShader.bind();
		mPositionFieldShader.uploadAttrib("o", mRes.O());
		//mPositionFieldShader.uploadAttrib("BC", mRes.BC_render);

		mLayers[Layers::BoundaryWireframe]->setChecked(true);
		mLayers[Layers::PositionField]->setChecked(true);
    });
//Rosy
    mSolveOrientationBtn = new Button(window, "Rosy", ENTYPO_ICON_FLASH);
    mSolveOrientationBtn->setBackgroundColor(Color(0, 0, 255, 25));
    mSolveOrientationBtn->setFlags(Button::Flags::ToggleButton);
    mSolveOrientationBtn->setChangeCallback([&](bool value) {
        mOptimizer->setOptimizeOrientations(value);
        mOptimizer->notify();

        mLayers[Layers::OrientationField]->setChecked(true);
        mLayers[Layers::PositionField]->setChecked(false);
    });	
//Posy
	mSolvePositionBtn = new Button(window, "Posy", ENTYPO_ICON_FLASH);
    mSolvePositionBtn->setBackgroundColor(Color(0, 0, 255, 25));
    mSolvePositionBtn->setFlags(Button::Flags::ToggleButton);
    mSolvePositionBtn->setChangeCallback([&](bool value) {

        mOptimizer->setOptimizePositions(value);
        mOptimizer->notify();

        mLayers[Layers::OrientationField]->setChecked(false);
        mLayers[Layers::PositionField]->setChecked(true);
    });
//Extraction
	mFinalEdgeTagging = new CheckBox(window, "Graph");
	mFinalEdgeTagging->setId("showedgetags");
	mFinalEdgeTagging->setChecked(false);

    mExtractBtn = new Button(window, "Extract", ENTYPO_ICON_FLASH);
    mExtractBtn->setBackgroundColor(Color(0, 255, 0, 25));
    mExtractBtn->setCallback([&]() {
		if(mRes.D3)
			mRes.edge_tagging3D();
		else
			mRes.edge_tagging2D();

		mRes.construct_Graph(mRes.m, mRes.g);
		//mRes.graph_extraction();
		mRes.graph_extraction_fast();
		mRes.composit_edges_colors(mRes.g, mRes.E_final_rend);

		mOrientationFieldShaderTri.bind();
		mOrientationFieldShaderTri.uploadAttrib("position", mRes.mV_final);
		mOrientationFieldShaderTri.uploadAttrib("q", mRes.mQ_final);
		mOrientationFieldShaderTri.uploadAttrib("n", mRes.mN_final);

		auto const &R = mRes.E_final_rend;
		mExtractionResultShader.bind();
		mExtractionResultShader.uploadAttrib("position", MatrixXf(R.block(0, 0, 3, R.cols())));
		mExtractionResultShader.uploadAttrib("color", MatrixXf(R.block(3, 0, 3, R.cols())));

		mFinalEdgeTagging->setChecked(true);
	});
	mExtractBtn_step = new Button(window, "Extract_step", ENTYPO_ICON_FLASH);
	mExtractBtn_step->setBackgroundColor(Color(0, 255, 0, 25));
	mExtractBtn_step->setCallback([&]() {
		mRes.graph_extraction();
		mRes.composit_edges_colors(mRes.g, mRes.E_final_rend);

		auto const &R = mRes.E_final_rend;
		mExtractionResultShader.bind();
		mExtractionResultShader.uploadAttrib("position", MatrixXf(R.block(0, 0, 3, R.cols())));
		mExtractionResultShader.uploadAttrib("color", MatrixXf(R.block(3, 0, 3, R.cols())));

		mFinalEdgeTagging->setChecked(true);
	});


 	//Config Layers
	PopupButton *openBtn3 = new PopupButton(window, "Config Layers");
	auto popup3 = openBtn3->popup();
	popup3->setLayout(new GroupLayout());
	Configlayers[Config_Layers::Alignment] = new CheckBox(popup3, "Boundary alignment");
	Configlayers[Config_Layers::Extrinsic] = new CheckBox(popup3, "Extrinsic smoothing");
	Configlayers[Config_Layers::Randomization] = new CheckBox(popup3, "Randomization");
	Configlayers[Config_Layers::Hierarchy] = new CheckBox(popup3, "Hierarchy");
	int ctr = 0;
	for (auto l : Configlayers) {
		l->setChecked(true);
		l->setId("configlayer." + std::to_string(ctr++));
	}
	Configlayers[Config_Layers::Extrinsic]->setChecked(true);
	//Render Layers
	PopupButton *openBtn = new PopupButton(window, "Render Layers");
	auto popup = openBtn->popup();
	popup->setLayout(new GroupLayout());

	mLayers[Layers::OrientationField] = new CheckBox(popup, "Orientation field");
	mLayers[Layers::PositionField] = new CheckBox(popup, "Position field");
	mLayers[Layers::BoundaryWireframe] = new CheckBox(popup, "Boundary wireframe");
	mLayers[Layers::Rosy_Out] = new CheckBox(popup, "Rosy_Out");

	ctr = 0;
	for (auto l : mLayers) {
		l->setChecked(false);
		l->setId("layer." + std::to_string(ctr++));
	}
//morphing

	PopupButton *MorphingBtn = new PopupButton(window, "Morphing");
	Popup *morphPopup = MorphingBtn->popup();
	morphPopup->setAnchorHeight(61);

	morphPopup->setLayout(new GroupLayout());

	mEdgeTagging = new CheckBox(morphPopup, "Coloring");
	mEdgeTagging->setId("showedgetags");
	mEdgeTagging->setChecked(false);
	mEdgeTagging->setCallback([&](bool value) {
		auto const &R = mRes.E_rend;
		mEdge_color_morphShader.bind();
		mEdge_color_morphShader.uploadAttrib("position", MatrixXf(R.block(0, 0, 3, R.cols())));
		mEdge_color_morphShader.uploadAttrib("color", MatrixXf(R.block(3, 0, 3, R.cols())));

		if(value){
			mLayers[Layers::PositionField]->setChecked(false);
			mLayers[Layers::BoundaryWireframe]->setChecked(false);
		}
	});
	Slider *slider2 = new Slider(morphPopup);
	slider2->setValue(0.0);
	auto cb = [&](Float value) {
		mRes.E_I_rend = (1 - value) * mRes.E_rend + value*mRes.E_O_rend;
	};
	cb(0.0f);
	slider2->setCallback(cb);
	slider2->setId("slider2");
//result morph
	mEdgeTaggingRanked = new CheckBox(morphPopup, "Ranked");
	mEdgeTaggingRanked->setId("showedgetags");
	mEdgeTaggingRanked->setChecked(false);
	mEdgeTaggingRanked->setCallback([&](bool value) {
		auto const &R = mRes.E_rend_o;
		mEdge_color_morphRankedShader.bind();
		mEdge_color_morphRankedShader.uploadAttrib("position", MatrixXf(R.block(0, 0, 3, R.cols())));
		mEdge_color_morphRankedShader.uploadAttrib("color", MatrixXf(R.block(3, 0, 3, R.cols())));

		if (value) {
			mLayers[Layers::PositionField]->setChecked(false);
			mLayers[Layers::BoundaryWireframe]->setChecked(false);
		}
	});
	Slider *slider3 = new Slider(morphPopup);
	slider3->setValue(0.0);
	auto cb2 = [&](Float value) {
		mRes.E_I_rend_o = (1 - value) * mRes.E_rend_o + value*mRes.E_O_rend_o;
	};
	cb2(0.0f);
	slider3->setCallback(cb2);
	slider3->setId("slider3");

//output
	mOutputBtn = new Button(window, "Output", ENTYPO_ICON_FLASH);
	mOutputBtn->setBackgroundColor(Color(0, 255, 0, 25));
	mOutputBtn->setCallback([&]() {
		char patho[300];
		sprintf(patho, "%s%s", filename.c_str(), "_graph.obj");
		write_Graph(patho, mRes.mV_final, mRes.mE_final);
		sprintf(patho, "%s%s", filename.c_str(), "_graph_opt.obj");
		write_Graph(patho, mRes.g);
	});
//layout
    performLayout();
}

Viewer::~Viewer() {
    mOptimizer->shutdown();
    delete mOptimizer;
}

bool Viewer::mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
                              int button, int modifiers) {
    if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
        if (mCamera.arcball.motion(p)) {
            //repaint();
        } else if (mTranslate) {
            Eigen::Matrix4f model, view, proj;
            computeCameraMatrices(model, view, proj);
            float zval = project(mRes.aabb().center().cast<float>(), view * model, proj, mSize).z();
            Eigen::Vector3f pos1 =
                unproject(Eigen::Vector3f(p.x(), mSize.y() - p.y(), zval),
                          view * model, proj, mSize);
            Eigen::Vector3f pos0 = unproject(
                Eigen::Vector3f(mTranslateStart.x(),
                                mSize.y() - mTranslateStart.y(), zval),
                view * model, proj, mSize);
            mCamera.modelTranslation =
                mCamera.modelTranslation_start + (pos1 - pos0);
            //repaint();
        }
    }
    return true;
}

bool Viewer::mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) {
    if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
        if (button == GLFW_MOUSE_BUTTON_1 && modifiers == 0) {
            mCamera.arcball.button(p, down);
        } else if (button == GLFW_MOUSE_BUTTON_2 ||
                   (button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_SHIFT)) {
            mCamera.modelTranslation_start = mCamera.modelTranslation;
            mTranslate = true;
            mTranslateStart = p;
        }
    }
    if (button == GLFW_MOUSE_BUTTON_1 && !down)
        mCamera.arcball.button(p, false);
    if (!down) {
        mTranslate = false;
    }
    return true;
}

bool Viewer::resizeEvent(const Vector2i &size) {
    mCamera.arcball.setSize(mSize);
    return true;
}

bool Viewer::scrollEvent(const Vector2i &p, const Eigen::Vector2f &rel) {
    if (!Screen::scrollEvent(p, rel)) {
        mCamera.zoom = std::max(0.1, mCamera.zoom * (rel.y() > 0 ? 1.1 : 0.9));
        //repaint();
    }
    return true;
}

void Viewer::drawContents() {
	glClearColor(.5, .5, .5, 1.0f);
	//glClearColor(1, 1, 1, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	mOptimizer->setAlignment(Configlayers[Config_Layers::Alignment]->checked());
	mOptimizer->setRandomization(Configlayers[Config_Layers::Randomization]->checked());
	mOptimizer->setHierarchy(Configlayers[Config_Layers::Hierarchy]->checked());
	mOptimizer->setExtrinsic(Configlayers[Config_Layers::Extrinsic]->checked());

	mRes.setScale(mScaleBox->value());

	if (!mOptimizer->active()) {
		if (mSolveOrientationBtn->pushed()) {
			mSolveOrientationBtn->setPushed(false);
		}
		if (mSolvePositionBtn->pushed()) {
			mSolvePositionBtn->setPushed(false);
		}
	}
	else if (!mOptimizer->hierarchy()) {
		if (mSolveOrientationBtn->pushed());
	}

	if (mSolveDatastructureBtn->pushed()) {
		mSolveDatastructureBtn->setPushed(false);
	}
	Eigen::Matrix4f model, view, proj;
	computeCameraMatrices(model, view, proj);
	Eigen::Matrix4f mvp = proj * view * model;
	Eigen::Vector4f civ =
		(view * model).inverse() * Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

	//if (mRes.D3) {
	//	mOrientationFieldShaderTet.bind();
	//	mOrientationFieldShaderTet.uploadAttrib("q", mRes.Q());
	//}

	//mPositionFieldShader.uploadAttrib("BC", mRes.BC_render);

	//glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LEQUAL);
	//glDisable(GL_BLEND);

	glPointSize(5);
	// This must be enabled, otherwise glLineWidth has no effect
	//glEnable(GL_LINE_SMOOTH);
	glLineWidth(1);
	if (mLayers[OrientationField]->checked()) {
		if (!mRes.D3)
		{
			mOrientationFieldShaderTri.bind();
			mOrientationFieldShaderTri.uploadAttrib("q", mRes.Q());
		}
		else
		{
			mOrientationFieldShaderTet.bind();
			mOrientationFieldShaderTet.uploadAttrib("q", mRes.Q());
		}


		auto &shader = mRes.D3 ? mOrientationFieldShaderTet : mOrientationFieldShaderTri;
		shader.bind();
		shader.setUniform("mvp", mvp);
		shader.setUniform("split", mSplit, false);
		//shader.setUniform("scale", mRes.ms.mAverageEdgeLength / 3);
		//shader.drawArray(GL_POINTS, 0, mRes.D3 ? mRes.tetCount() : mRes.faceCount());
		shader.drawArray(GL_POINTS, 0, mRes.D3 ? mRes.vertexCount() : mRes.vertexCount());
	}

	if (mLayers[PositionField]->checked()) {
		mPositionFieldShader.bind();
		mPositionFieldShader.uploadAttrib("o", mRes.O());
		mPositionFieldShader.bind();
		mPositionFieldShader.setUniform("mvp", mvp);
		mPositionFieldShader.setUniform("split", mSplit, false);
		mPositionFieldShader.drawArray(GL_POINTS, 0, mRes.vertexCount());
	}

	if (mLayers[BoundaryWireframe]->checked()) {
		auto &shader = mGraphShader;
		shader.bind();
		shader.setUniform("split", mSplit, false);
		shader.setUniform("mvp", mvp);
		shader.drawArray(GL_LINES, 0, mRes.E_input_rend.cols());
	}

	if (mEdgeTagging->checked()) {
		mEdge_color_morphShader.bind();
		mEdge_color_morphShader.uploadAttrib("position", MatrixXf(mRes.E_I_rend.block(0, 0, 3, mRes.E_I_rend.cols())));
		mEdge_color_morphShader.uploadAttrib("color", MatrixXf(mRes.E_I_rend.block(3, 0, 3, mRes.E_I_rend.cols())));
		auto &shader = mEdge_color_morphShader;
		shader.bind();
		shader.setUniform("split", mSplit, false);
		shader.setUniform("mvp", mvp);
		shader.drawArray(GL_LINES, 0, mRes.E_I_rend.cols());
	}
	if (mEdgeTaggingRanked->checked()) {
		mEdge_color_morphRankedShader.bind();
		mEdge_color_morphRankedShader.uploadAttrib("position", MatrixXf(mRes.E_I_rend_o.block(0, 0, 3, mRes.E_I_rend_o.cols())));
		mEdge_color_morphRankedShader.uploadAttrib("color", MatrixXf(mRes.E_I_rend_o.block(3, 0, 3, mRes.E_I_rend_o.cols())));
		auto &shader = mEdge_color_morphRankedShader;
		shader.bind();
		shader.setUniform("split", mSplit, false);
		shader.setUniform("mvp", mvp);
		shader.drawArray(GL_LINES, 0, mRes.E_I_rend_o.cols());
	}

	 if (mFinalEdgeTagging->checked())
	 {
	 	auto &shader = mExtractionResultShader;
	 	shader.bind();
	 	//shader.setUniform("split", mSplit, false);
	 	shader.setUniform("mvp", mvp);
	 	shader.drawArray(GL_LINES, 0, mRes.E_final_rend.cols());
	 }
	 if (mLayers[Rosy_Out]->checked()) {
		 mOrientationFieldShaderTri.bind();
		 mOrientationFieldShaderTri.uploadAttrib("q", mRes.mQ_final);

		 auto &shader = mOrientationFieldShaderTri;
		 shader.bind();
		 shader.setUniform("mvp", mvp);
		 shader.setUniform("split", mSplit, false);
		 shader.drawArray(GL_POINTS, 0, mRes.mV_final.cols());
	 }

}

bool Viewer::keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (Screen::keyboardEvent(key, scancode, action, modifiers))
        return true;
    if (action != GLFW_PRESS)
        return false;
     return false;
}

void Viewer::computeCameraMatrices(Eigen::Matrix4f &model,
                                   Eigen::Matrix4f &view,
                                   Eigen::Matrix4f &proj) {
    view = lookAt(mCamera.eye, mCamera.center, mCamera.up);

    float fH = std::tan(mCamera.viewAngle / 360.0f * M_PI) * mCamera.dnear;
    float fW = fH * (float) mSize.x() / (float) mSize.y();

    proj = frustum(-fW, fW, -fH, fH, mCamera.dnear, mCamera.dfar);
    model = mCamera.arcball.matrix();

	model = model * scale(Eigen::Vector3f::Constant(mCamera.zoom * mCamera.modelZoom));
	model = model * translate(mCamera.modelTranslation);
}


