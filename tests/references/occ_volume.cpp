// Independent adaptive OCCT mass properties of a BREP imported by FreeCAD.
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Shape.hxx>
#include <cmath>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    BRep_Builder builder;
    TopoDS_Shape shape;
    if (!BRepTools::Read(shape, argv[1], builder) || shape.IsNull()) return 3;
    std::cout << std::setprecision(17) << "[";
    for (const double eps : {1.e-9, 1.e-12}) {
        GProp_GProps props;
        const double error = BRepGProp::VolumePropertiesGK(shape, props, eps, true, true);
        if (!std::isfinite(error) || error < 0 || !std::isfinite(props.Mass())) return 4;
        if (eps == 1.e-12) std::cout << ",";
        std::cout << "{\"tolerance\":" << eps << ",\"volume_mm3\":" << props.Mass()
                  << ",\"estimated_error\":" << error << "}";
    }
    std::cout << "]\n";
}
