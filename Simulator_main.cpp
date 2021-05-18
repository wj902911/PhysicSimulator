#include <core/shared/helper.h>
//#include <vis/cloth2/ClothVisualizer.h>
//#include <vis/goo2/GooVisualizer.h>
#include <vis/fluid/FluidVisualizer.h>

int main()
{
	//Eigen::MatrixX3d V;
	//Eigen::MatrixX3i F;
	//Eigen::MatrixX3d V2;
	//Eigen::MatrixX3i F2;
	//std::tie(V, F) = loadOBJ("../assets/cloth/rect-coarse.obj");
	//std::tie(V2, F2) = loadOBJ("../assets/cloth/longhorn2.obj");

	//cloth2::ClothVisualizer vis = cloth2::ClothVisualizer();
	//goo2::GooVisualizer vis = goo2::GooVisualizer();
	fluid::FluidVisualizer vis = fluid::FluidVisualizer();
	//vis.attachMesh(V, F);
	//vis.attachBody(V2, F2, 0.05);
	vis.run();
}