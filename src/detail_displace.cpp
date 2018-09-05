#include "smooth.h"

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <thread>
#include <exception>
#include "cmdline.h"
#include "timer.h"

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;



cmdline::parser cmd_parser;
void cmd_parser_build()
{
	cmd_parser.add<std::uint32_t>("f",  0, "Number of filter iterations", false);
	cmd_parser.add<std::string>("d",  0, "Detail mesh file", true);
	cmd_parser.add<std::string>("t",  0, "Target mesh file", false);
	cmd_parser.add<std::string>("df", 0, "Filtered detail mesh file", true);
	cmd_parser.add<std::string>("tf", 0, "Filtered target mesh file", true);
	cmd_parser.add<std::string>("o",  0, "Output mesh result filename", true, "output.obj");
}


int main(int argc, char* argv[])
{
	//
	// Parse command line
	//
	cmd_parser_build();
	cmd_parser.parse_check(argc, argv);


	//
	// Set variables
	//
	const std::string& detail_mesh_filename 		= cmd_parser.get<std::string>("d");
	const std::string& detail_base_mesh_filename	= cmd_parser.get<std::string>("df");
	const std::string& target_mesh_filename 		= cmd_parser.get<std::string>("t");
	const std::string& target_base_mesh_filename 	= cmd_parser.get<std::string>("tf");
	const std::string& output_mesh_filename 		= cmd_parser.get<std::string>("o");
	const uint32_t number_of_filter_iterations		= cmd_parser.get<uint32_t>("f");


	Smooth::MyMesh detail_mesh;	
	//Smooth::MyMesh target_mesh;
	Smooth::MyMesh detail_base_mesh;
	Smooth::MyMesh target_base_mesh;

	timer t;
	t.start();

	std::vector<std::thread> thread;
	try
	{
		std::cout << "Loading: " << detail_mesh_filename << '\n';
		thread.push_back(std::thread([&] { OpenMesh::IO::read_mesh(detail_mesh, detail_mesh_filename); }));

		//std::cout << "Loading: " << target_mesh_filename << '\n';
		//thread.push_back(std::thread([&]() { OpenMesh::IO::read_mesh(target_mesh, target_mesh_filename); }));

		std::cout << "Loading: " << detail_base_mesh_filename << '\n';
		thread.push_back(std::thread([&] { OpenMesh::IO::read_mesh(detail_base_mesh, detail_base_mesh_filename); }));

		std::cout << "Loading: " << target_base_mesh_filename << '\n';
		thread.push_back(std::thread([&] { OpenMesh::IO::read_mesh(target_base_mesh, target_base_mesh_filename); }));
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

	for (auto& t : thread) 
	 	t.join();


	t.stop();
	t.print_interval_msec("Load meshes time: ");

	if (detail_mesh.n_vertices() != target_base_mesh.n_vertices() || detail_base_mesh.n_vertices() != target_base_mesh.n_vertices())
	{
		std::cerr << "Error: Vertex count does not match. The meshes must have the same topology" << std::endl;
		return EXIT_FAILURE;
	}

	Smooth::MyMesh::VertexIter v_it_d(detail_mesh.vertices_begin()), v_end_d(detail_mesh.vertices_end());
	Smooth::MyMesh::VertexIter v_it_db(detail_base_mesh.vertices_begin()), v_end_db(detail_base_mesh.vertices_end());
	Smooth::MyMesh::VertexIter v_it_tb(target_base_mesh.vertices_begin()), v_end_tb(target_base_mesh.vertices_end());

	for (; v_it_tb != v_end_tb; ++v_it_d, ++v_it_db, ++v_it_tb) 
	{
		const auto vec3 = detail_mesh.point(*v_it_d) - detail_base_mesh.point(*v_it_db) + target_base_mesh.point(*v_it_tb);
		detail_mesh.set_point(*v_it_d, vec3);
	}

	
	std::cout << "Saving : " << output_mesh_filename << '\n';
	if (!OpenMesh::IO::write_mesh(detail_mesh, output_mesh_filename))
	{
		std::cerr << "Error: cannot write mesh to " << output_mesh_filename << std::endl;
		return false;
	}

	return EXIT_SUCCESS;
}
