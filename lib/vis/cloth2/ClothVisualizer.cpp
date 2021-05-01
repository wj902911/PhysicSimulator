#include "ClothVisualizer.h"
#include "ClothHook.h"
#include <core/cloth2/ClothCore.h>
#include <core/cloth2/SimParameters.h>
#include <igl/unproject_onto_mesh.h>

namespace cloth2 {

const char* kMouseModes[] = { "Rotate Camera", "Drag Cloth", "Air Jet", "Saber", "Bevo Launcher" };

ClothVisualizer::ClothVisualizer()
    :mouseMode_(MouseMode::MM_ROTATECAMERA)
{
    core_.reset(new ClothCore);
    hook_.reset(new ClothHook(core_.get()));
    init(core_.get(), hook_.get());
}

ClothVisualizer::~ClothVisualizer()
{
    // Must kill simulation thread before deleting core_
    hook_.reset(); // called killSimThread in the destructor
    core_.reset(); // No necessary but make it explicit.
}

void ClothVisualizer::setupViewer()
{
    viewer_.data().set_face_based(true);
    viewer_.core().is_animating = true;

    viewer_.callback_mouse_up = [=](Viewer& viewer, int button, int mod) {
        return this->mouseReleased(viewer, button, mod);
    };
    viewer_.callback_mouse_move = [=](Viewer& viewer, int mouse_x, int mouse_y) {
        return this->mouseMoved(viewer, mouse_x, mouse_y);
    };
}

bool ClothVisualizer::mouseCallback(Viewer& viewer, int button, int modifier)
{
    mouse_pressed_ = true;
    if (button != static_cast<int>(Viewer::MouseButton::Left))
        return false;
    if (mouseMode_ == MM_ROTATECAMERA)
        return false;

    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    int fid;
    Eigen::Vector3f bc;
    MouseEvent me;
    me.mode = mouseMode_;
    bool ret = true;

    hook_->lockRenderer(); // Protect the result of hook_->getV()
    if (mouseMode_ == MouseMode::MM_AIRJET) {
        handleJetEvent(viewer, me);
    } else if (mouseMode_ == MouseMode::MM_DRAGCLOTH) {
        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y),
                                     viewer.core().view,
                                     viewer.core().proj,
                                     viewer.core().viewport,
                                     hook_->getV(), hook_->getF(),
                                     fid, bc)) {
            int bestvert = -1;
            double bestcoord = 2.0;
            for (int j = 0; j < 3; j++) {
                if (bc[j] < bestcoord) {
                    bestcoord = bc[j];
                    bestvert = j;
                }
            }
            me.type = MouseEvent::ME_CLICKED;
            me.vertex = hook_->getF()(fid, bestvert);

            Eigen::Vector3f proj;
            Eigen::Vector3f pt = hook_->getV().row(me.vertex).cast<float>();
            Eigen::Matrix4f modelview = viewer.core().view;
            proj = igl::project(pt, modelview, viewer.core().proj, viewer.core().viewport);

            clickedz_ = proj[2];
            Eigen::Vector3f pos;
            pos << x, y, clickedz_;
            Eigen::Vector3f unproj = igl::unproject(pos,
                                                    modelview,
                                                    viewer.core().proj,
                                                    viewer.core().viewport);
            me.pos = unproj.cast<double>();
            ret = true;
        } else {
            me.type = MouseEvent::ME_RELEASED;
            ret = false;
        }
    } else if (mouseMode_ == MouseMode::MM_SABER) {
        handleSaberEvent(viewer, me);
    } else if (mouseMode_ == MouseMode::MM_BEVO) {
        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
        Eigen::Vector3f pos;
        pos << x, y, 0.5;
        Eigen::Vector3f wrdpos;
        igl::unproject(pos, viewer.core().view, viewer.core().proj, viewer.core().viewport, wrdpos);

        // me.pos = viewer.core().camera_eye.cast<double>();
        me.pos = viewer.core().view.cast<double>().inverse().col(3).segment<3>(0);
        me.dir = (wrdpos.cast<double>() - me.pos).normalized();
        me.type = MouseEvent::ME_CLICKED;
    }
    hook_->unlockRenderer();

    hook_->addMouseEvent(me);
    return ret;
}

bool ClothVisualizer::mouseReleased(Viewer& viewer, int button, int mod)
{
    mouse_pressed_ = false;
    if (mouseMode_ == MM_ROTATECAMERA)
        return false;
    MouseEvent me;
    me.type = MouseEvent::ME_RELEASED;
    hook_->addMouseEvent(me);
    return false;
}

bool ClothVisualizer::mouseMoved(Viewer& viewer, int mouse_x, int mouse_y)
{
    if (mouseMode_ == MM_ROTATECAMERA)
        return false;
    if (!mouse_pressed_)
        return false;

    MouseEvent me;
    me.mode = mouseMode_;
    if (mouseMode_ == MouseMode::MM_AIRJET) {
        handleJetEvent(viewer, me);
    } else if (mouseMode_ == MouseMode::MM_DRAGCLOTH) {
        me.type = MouseEvent::ME_DRAGGED;
        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
        Eigen::Vector3d pos(x, y, clickedz_);
        igl::unproject(pos, viewer.core().view, viewer.core().proj, viewer.core().viewport, me.pos);
    } else if (mouseMode_ == MouseMode::MM_SABER) {
        handleSaberEvent(viewer, me);
    } else if (mouseMode_ == MouseMode::MM_BEVO) {
        return false;
    }
    hook_->addMouseEvent(me);
    return false;
}

bool ClothVisualizer::drawGUI()
{
    auto param = core_->getPointerToSimParameters();
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputFloat("Timestep", &param->dt, 0, 0, "%.3f");
        ImGui::InputInt("Constraint Iters", &param->constraintIters);
    }
    ImGui::Combo("Mouse Mode", (int*)&mouseMode_, kMouseModes, IM_ARRAYSIZE(kMouseModes));
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Gravity Enabled", &param->gravityEnabled);
        ImGui::InputFloat("Gravity G", &param->gravityG, 0, 0, "%.3f");
        ImGui::Checkbox("Pins Enabled", &param->pinEnabled);
        ImGui::InputFloat("Pin Weight", &param->pinWeight, 0, 0, "%.3f");
        ImGui::Checkbox("Stretching Enabled", &param->stretchEnabled);
        ImGui::InputFloat("Stretching Weight", &param->stretchWeight, 0, 0, "%.3f");
        ImGui::Checkbox("Bending Enabled", &param->bendingEnabled);
        ImGui::InputFloat("Bending Weight", &param->bendingWeight, 0, 0, "%.3f");
        ImGui::Checkbox("Collisions Enabled", &param->collisionsEnabled);
        ImGui::InputFloat("Drag Weight", &param->dragWeight, 0, 0, "%.3f");
        ImGui::InputFloat("Jet Strength", &jetStrength_, 0, 0, "%.3f");
        ImGui::InputFloat("Jet Falloff", &jetAtten_, 0, 0, "%.3f");
        ImGui::InputFloat("Launch Speed", &param->launchSpeed, 0, 0, "%.3f");
        ImGui::InputFloat("Launch Spinning Speed", &param->launchSpinningSpeed, 0, 0, "%.3f");
    }
    if(ImGui::Button("Deploy Bevo"))
    {
        param->shouldDeployBevo = true;
    }

    return IglVisualizer::drawGUI();
}

void ClothVisualizer::attachMesh(Eigen::Ref<Eigen::MatrixX3d> V,
                                 Eigen::Ref<Eigen::MatrixX3i> F,
                                 double scale)
{
    core_->attachMesh(V, F, scale);
    resetSimulation();
}

void ClothVisualizer::attachBody(Eigen::Ref<Eigen::MatrixX3d> V,
                                 Eigen::Ref<Eigen::MatrixX3i> F,
                                 double scale)
{
    core_->attachBody(V, F, scale);
    resetSimulation();
}

bool ClothVisualizer::unproject_onto_mesh(Viewer& viewer,
                                          Eigen::Vector3f& scrpos,
                                          Eigen::Vector3f& wrdpos,
                                          int& fid, Eigen::Vector3f& bc) const
{
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    const auto& F = hook_->getF();
    const auto& V = hook_->getV();
    Eigen::Matrix4f modelview = viewer.core().view;

    if (!igl::unproject_onto_mesh(Eigen::Vector2f(x, y),
                                  modelview,
                                  viewer.core().proj,
                                  viewer.core().viewport,
                                  V, F,
                                  fid, bc))
        return false;
    Eigen::Vector3i f = F.row(fid);
    Eigen::Vector3f pt = (V.row(f(0)) * bc(0) + V.row(f(1)) * bc(1) + V.row(f(2)) * bc(2)).cast<float>();
    Eigen::Vector3f proj = igl::project(pt, modelview, viewer.core().proj, viewer.core().viewport);

    Eigen::Vector3f pos;
    pos << x, y, proj(2);
    igl::unproject(pos, viewer.core().view, viewer.core().proj, viewer.core().viewport, wrdpos);
    return true;
}

bool ClothVisualizer::handleJetEvent(Viewer& viewer, MouseEvent& me) const
{
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    Eigen::Vector3f pos;
    pos << x, y, 0.5;
    Eigen::Vector3f wrdpos;
    igl::unproject(pos, viewer.core().view, viewer.core().proj, viewer.core().viewport, wrdpos);

    // me.pos = viewer.core().camera_eye.cast<double>();
    me.pos = viewer.core().view.cast<double>().inverse().col(3).segment<3>(0);
    me.dir = (wrdpos.cast<double>() - me.pos).normalized();
    me.jetStrength = jetStrength_;
    me.jetAtten = jetAtten_;
    return true;
}

bool ClothVisualizer::handleSaberEvent(Viewer& viewer, MouseEvent& me) const
{
    Eigen::Vector3f scrpos;
    Eigen::Vector3f wrdpos;
    int fid;
    Eigen::Vector3f bc;
    if (unproject_onto_mesh(viewer, scrpos, wrdpos, fid, bc)) {
        me.pos = viewer.core().view.cast<double>().inverse().col(3).segment<3>(0);
        me.dir = (wrdpos.cast<double>() - me.pos).normalized().cast<double>();
    }
    return true;
}

}
