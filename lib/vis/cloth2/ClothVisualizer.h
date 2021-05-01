#include "../IglVisualizer.h"
#include <memory>
#include <thread>

namespace cloth2 {

class ClothCore;
class ClothHook;
enum MouseMode : int;
struct MouseEvent;

class ClothVisualizer : public IglVisualizer {
protected:
    std::shared_ptr<ClothHook> hook_;
    std::shared_ptr<ClothCore> core_;

public:
    ClothVisualizer();
    ~ClothVisualizer();

    virtual void setupViewer() override;

    virtual bool mouseCallback(Viewer& viewer, int button, int mod) override;

    virtual bool mouseReleased(Viewer& viewer, int button, int mod);
    virtual bool mouseMoved(Viewer& viewer, int mouse_x, int mouse_y);

    virtual bool drawGUI() override;

    std::shared_ptr<ClothCore> getCore() { return core_; }

    void attachMesh(Eigen::Ref<Eigen::MatrixX3d> V,
                    Eigen::Ref<Eigen::MatrixX3i> F,
                    double scale = 50.0);

    void attachBody(Eigen::Ref<Eigen::MatrixX3d> V,
                    Eigen::Ref<Eigen::MatrixX3i> F, double scale = 1.0);

protected:
    double clickedz_;
    MouseMode mouseMode_;
    bool mouse_pressed_ = true;
    float jetStrength_ = 1.0;
    float jetAtten_ = 0.0;

    bool unproject_onto_mesh(Viewer& viewer,
                             Eigen::Vector3f& scrpos,
                             Eigen::Vector3f& wrdpos,
                             int& fid, Eigen::Vector3f& bc) const;

    bool handleSaberEvent(Viewer& viewer, MouseEvent& me) const;
    bool handleJetEvent(Viewer& viewer, MouseEvent& me) const;
};

}
