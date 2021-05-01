#include "ClothHook.h"
#include <core/cloth2/ClothCore.h>
#include <core/cloth2/SimParameters.h>
#include <core/cloth2/SceneObjects.h>
#include <core/shared/VectorMath.h>
#include <igl/file_dialog_open.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

using namespace Eigen;

namespace cloth2 {

ClothHook::ClothHook(ClothCore* core)
    : PhysicsHook(core)
    , core_(core)
    , params_(*core->getPointerToSimParameters())
{
}

ClothHook::~ClothHook()
{
}

void ClothHook::tick()
{
    mouse_mutex_.lock();
    std::vector<int> face_to_cut;
    for (MouseEvent me : mouse_events_) {
        params_.dragEnabled = false;
        if (me.mode == MouseMode::MM_DRAGCLOTH) {
            if (me.type == MouseEvent::ME_CLICKED) {
                params_.dragEnabled = true;
                core_->hold(me.vertex, me.pos);
            }
            if (me.type == MouseEvent::ME_RELEASED) {
                core_->releaseHold();
            }
            if (me.type == MouseEvent::ME_DRAGGED) {
                params_.dragEnabled = true;
                core_->updateHold(me.pos);
            }
        } else if (me.mode == MouseMode::MM_AIRJET) {
            auto jet = core_->queryJetSource(jet_id_);
            if (!jet && me.type != MouseEvent::ME_RELEASED) {
                // Ensure the jet source
                jet_id_ = core_->addJetSource(me.pos, me.dir,
                                              me.jetStrength, me.jetAtten);
                jet = core_->queryJetSource(jet_id_);
            }
            if (jet) {
                jet->updateParameter(me.jetStrength, me.jetAtten);
                jet->updateConfiguration(me.pos, me.dir);
                if (me.type == MouseEvent::ME_RELEASED) {
                    jet->disable();
                } else {
                    jet->enable();
                }
            }
        } else if (me.mode == MouseMode::MM_SABER) {
            std::vector<igl::Hit> hits;
            igl::ray_mesh_intersect(me.pos, me.dir,
                                    getV(), getF(),
                                    hits);
            std::vector<int> face_hits;
            face_hits.reserve(hits.size());
            for (const auto& hit : hits)
                face_hits.emplace_back(hit.id);
            face_to_cut.insert(face_to_cut.end(), face_hits.begin(), face_hits.end());
        } else if (me.mode == MouseMode::MM_BEVO) {
            if (me.type == MouseEvent::ME_CLICKED) {
                auto orientation = VectorMath::randomPointOnSphere();
                auto omega = VectorMath::randomPointOnSphere();
                core_->launchBevo(me.pos, me.dir, orientation, omega);
            }
        }
    }
    core_->cutFaces(face_to_cut);
    mouse_events_.clear();
    mouse_mutex_.unlock();
    if (params_.shouldDeployBevo)
        core_->deployBevo();
}

void ClothHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu)
{
}

void ClothHook::addMouseEvent(const MouseEvent& me)
{
    mouse_mutex_.lock();
    mouse_events_.emplace_back(me);
    mouse_mutex_.unlock();
}

}
