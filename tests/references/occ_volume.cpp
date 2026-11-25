// Independent adaptive OCCT mass properties of a BREP imported by FreeCAD.
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Shape.hxx>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

int main(int argc, char** argv) {
    const bool centroid = argc == 3 && std::string_view(argv[2]) == "--centroid";
    if (argc != 2 && !centroid) return 2;
    BRep_Builder builder;
    TopoDS_Shape shape;
    if (!BRepTools::Read(shape, argv[1], builder) || shape.IsNull()) return 3;
    std::cout << std::setprecision(17) << "[";
    for (const double eps : {1.e-9, 1.e-12}) {
        GProp_GProps props;
        const double error = BRepGProp::VolumePropertiesGK(shape, props, eps, true, true, centroid);
        if (!std::isfinite(error) || error < 0 || !std::isfinite(props.Mass())) return 4;
        GProp_GProps surface;
        const double area_error = BRepGProp::SurfaceProperties(shape, surface, eps);
        if (!std::isfinite(area_error) || area_error < 0 || !std::isfinite(surface.Mass())) return 5;
        if (eps == 1.e-12) std::cout << ",";
        std::cout << "{\"tolerance\":" << eps << ",\"volume_mm3\":" << props.Mass()
                  << ",\"estimated_error\":" << error << ",\"surface_mm2\":" << surface.Mass()
                  << ",\"surface_estimated_error\":" << area_error;
        if (centroid) {
            const auto c = props.CentreOfMass();
            std::cout << ",\"centroid_mm\":[" << c.X() << "," << c.Y() << "," << c.Z() << "]";
        }
        std::cout << "}";
    }
    std::cout << "]\n";
}
