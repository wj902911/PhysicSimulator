#include "../PhysicsHook.h"
#include <deque>
#include <igl/opengl/glfw/Viewer.h>

namespace cloth2 {

struct SimParameters;
class ClothCore;

enum MouseMode : int
{
    MM_ROTATECAMERA = 0,
    MM_DRAGCLOTH = 1,
    MM_AIRJET = 2,
    MM_SABER = 3,
    MM_BEVO = 4,
};


struct MouseEvent {
    enum METype {
        ME_CLICKED,
        ME_RELEASED,
        ME_DRAGGED
    };

    METype type;
    MouseMode mode;
    int vertex;
    Eigen::Vector3d pos;
    Eigen::Vector3d dir;
    
    double jetAtten;
    double jetStrength;
};

class ClothHook : public PhysicsHook {
public:
    ClothHook(ClothCore* core);
    ~ClothHook();

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu) override;
    virtual void tick() override;

    void lockRenderer() { render_mutex.lock(); }
    void unlockRenderer() { render_mutex.unlock(); }

    const Eigen::MatrixXd& getV() const { return renderQ_; }
    const Eigen::MatrixXi& getF() const { return renderF_; }

    void addMouseEvent(const MouseEvent& me);

private:
    ClothCore* core_;
    SimParameters& params_;

    bool launch_ = false;
    Eigen::Vector3d launch_pos_;
    Eigen::Vector3d launch_dir_;

    std::mutex mouse_mutex_;
    std::vector<MouseEvent> mouse_events_;
    
    uint32_t jet_id_ = 0;
};

}
