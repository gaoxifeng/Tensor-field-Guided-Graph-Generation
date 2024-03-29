#pragma once

#include "common.h"
#include "optimizer.h"
#include <nanogui/nanogui.h>
#include <nanogui/glutil.h>

using nanogui::Alignment;
using nanogui::GLFramebuffer;
using nanogui::Screen;
using nanogui::Window;
using nanogui::Label;
using nanogui::Button;
using nanogui::frustum;
using nanogui::lookAt;
using nanogui::project;
using nanogui::scale;
using nanogui::translate;
using nanogui::unproject;
using nanogui::Orientation;
using nanogui::Arcball;
using nanogui::GLShader;
using nanogui::BoxLayout;
using nanogui::GroupLayout;
using nanogui::Slider;
using nanogui::PopupButton;
using nanogui::Popup;
using nanogui::Color;
using nanogui::CheckBox;
using nanogui::TextBox;
using nanogui::IntBox;
using nanogui::FloatBox;
using nanogui::ImagePanel;
using nanogui::VScrollPanel; 

class Viewer : public Screen {
public:
    Viewer(std::string &filename, bool fullscreen);
    ~Viewer();

protected:
    void drawContents();

    void computeCameraMatrices(Eigen::Matrix4f &model, Eigen::Matrix4f &view,
                               Eigen::Matrix4f &proj);

    bool resizeEvent(const Vector2i &size);

    bool mouseButtonEvent(const Vector2i &p, int button, bool down,
                          int modifiers);

    bool mouseMotionEvent(const Vector2i &p, const Vector2i &rel, int button,
                          int modifiers);

    bool keyboardEvent(int key, int scancode, int action, int modifiers);

    bool scrollEvent(const Vector2i &p, const Eigen::Vector2f &rel);

protected:
    struct CameraParameters {
        Arcball arcball;
        float zoom = 1.0f, viewAngle = 45.0f;
        float dnear = 0.05f, dfar = 100.0f;
        Eigen::Vector3f eye = Eigen::Vector3f(0.0f, 0.0f, 5.0f);
        Eigen::Vector3f center = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
        Eigen::Vector3f up = Eigen::Vector3f(0.0f, 1.0f, 5.0f);
        Eigen::Vector3f modelTranslation = Eigen::Vector3f::Zero();
        Eigen::Vector3f modelTranslation_start = Eigen::Vector3f::Zero();
        float modelZoom = 1.0f;
    };

	std::vector<std::pair<int, std::string>> mExampleImages;
	std::string mFilename;

	enum Config_Layers {
		Alignment = 0,
		Extrinsic,
		Randomization,
		Hierarchy,
		Config_LayerCount
	};
    enum Layers {
        OrientationField = 0,
        PositionField,
        BoundaryWireframe,
		Rosy_Out,
        LayerCount
    };
	enum Extraction_Condition {
		ReColor=0,
		Quadric,
		Splitting,
		Triangles,
		Doublets,
		Decompose,
		ExtractLayerCount
	};

    /* Mesh data */
    MultiResolutionHierarchy mRes;
    Optimizer *mOptimizer;

    /* UI */
	Button *mSolveDatastructureBtn;
    Button *mSolveOrientationBtn;
    Button *mSolvePositionBtn;
	CheckBox *Configlayers[Config_LayerCount];
    Button *mExtractBtn;
	Button *mExtractBtn_step;
	CheckBox *mFinalEdgeTagging;
    CheckBox *mLayers[LayerCount];
	FloatBox<Float> *mScaleBox;

	CheckBox *mEdgeTagging;
	CheckBox *mEdgeTaggingRanked;

	Button *mOutputBtn;

    /* Visualization */
    GLShader mOrientationFieldShaderTri, mOrientationFieldShaderTet;
    GLShader mPositionFieldShader;

    GLShader mGraphShader;
	GLShader mEdge_color_morphShader;
	GLShader mEdge_color_morphRankedShader;
    GLShader mExtractionResultShader;


    Vector4f mBaseColor, mSpecularColor;
    Vector4f mBaseColorBoundary, mSpecularColorBoundary;
    CameraParameters mCamera;
    Vector2i mTranslateStart;
    Vector4f mSplit;
    Vector3f mLightPosition;
    bool mTranslate;
};
