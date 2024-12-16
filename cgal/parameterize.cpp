#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

// border methods
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>

// parameterization methods
#include <CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>       Kernel;
typedef Kernel::Point_2                      Point_2;
typedef Kernel::Point_3                      Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>  SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor      face_descriptor;

namespace SMP = CGAL::Surface_mesh_parameterization;



int main(int argc, char** argv)
{
    // read mesh from stdin
    SurfaceMesh sm;
    std::cin >> sm;
    int nv = sm.number_of_vertices();
    if (nv == 0) {
        std::cerr << "Error: no vertex.\n";
        return EXIT_FAILURE;
    }

    // choose method
    char method = 'M';
    double param1 = 0;
    if (argc > 1)
        method = argv[1][0];
    if (argc > 2)
    	param1 = atof(argv[2]);


    // a halfedge on the border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

    // the 2D points of the uv parametrisation will be written into this map
    typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>  UV_pmap;
    UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("v:uv").first;

    // border parametrizers
    typedef SMP::Circular_border_arc_length_parameterizer_3<SurfaceMesh>  Border_A;
    typedef SMP::Two_vertices_parameterizer_3<SurfaceMesh>  Border_F;

    // surface parametrizers
    typedef SMP::Discrete_authalic_parameterizer_3<SurfaceMesh, Border_A> Parameterizer_AA;
    typedef SMP::Barycentric_mapping_parameterizer_3<SurfaceMesh, Border_A> Parameterizer_BA;
    typedef SMP::Discrete_conformal_map_parameterizer_3<SurfaceMesh, Border_A> Parameterizer_CA;
    typedef SMP::Mean_value_coordinates_parameterizer_3<SurfaceMesh, Border_A> Parameterizer_MA;
    typedef SMP::LSCM_parameterizer_3<SurfaceMesh, Border_F> Parameterizer_LF;
    typedef SMP::ARAP_parameterizer_3<SurfaceMesh, Border_F> Parameterizer_RF;

    // run parameterize
    SMP::Error_code err;
    if (method == 'A')
        err = SMP::parameterize(sm, Parameterizer_AA(), bhd, uv_map);
    else if (method == 'B')
        err = SMP::parameterize(sm, Parameterizer_BA(), bhd, uv_map);
    else if (method == 'C')
        err = SMP::parameterize(sm, Parameterizer_CA(), bhd, uv_map);
    else if (method == 'M')
        err = SMP::parameterize(sm, Parameterizer_MA(), bhd, uv_map);
    else if (method == 'L')
        err = SMP::parameterize(sm, Parameterizer_LF(), bhd, uv_map);
    else if (method == 'R')
        err = SMP::parameterize(sm, Parameterizer_RF(param1), bhd, uv_map);

    // error handling
    if (err != SMP::OK) {
        std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
    }

    // output parameterization to stdout
    if (method == 'R')
    	std::cout << "-----";

    BOOST_FOREACH(vertex_descriptor vd, sm.vertices()){
    	std::cout << get(uv_map, vd) << std::endl;
    }
    // output file
    //std::ofstream out("result.off");
    //SMP::IO::output_uvmap_to_off(sm, bhd, uv_map, out);
    return EXIT_SUCCESS;
}
