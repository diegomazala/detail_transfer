#include "smooth.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;


bool save(Smooth::MyMesh& mesh, int it, const std::string& filename)
{
	fs::path output_filename(filename);
	std::stringstream out;
	out << output_filename.parent_path() << '/' << output_filename.stem() << '_' << it << output_filename.extension();
	if (!OpenMesh::IO::write_mesh(mesh, out.str()))
	{
		std::cerr << "Error: cannot write mesh to " << out.str() << std::endl;
		return false;
	}
	else
	{
		std::cerr << "Ok   : Saved mesh to " << out.str() << std::endl;
		return true;
	}
}


int main(int argc, char **argv) 
{
	//return Smooth::main_test(argc, argv);

	Smooth::MyMesh mesh;
    
	Smooth::Mode mode = static_cast<Smooth::Mode>(atoi(argv[1]));
	uint32_t numberOfIterations = atoi(argv[2]);
	const std::string inputFilename = argv[3];
	const std::string outputFilename = argv[4];

	switch(mode)
	{
		case Smooth::Mode::Laplacian   			: std::cout << "Smooth::Mode::Laplacian\n";    			break;
		case Smooth::Mode::Taubin   			: std::cout << "Smooth::Mode::Taubin\n";    			break;
		case Smooth::Mode::LaplacianCotanWeight : std::cout << "Smooth::Mode::LaplacianCotanWeight\n";  break;
		default    								: std::cout << "Unknown option. Abort\n";  				return 1;
	}

	std::cout << "Loading: " << inputFilename << '\n';
	if (!OpenMesh::IO::read_mesh(mesh, inputFilename))
	{
		std::cerr << "Error: Cannot read mesh from " << inputFilename << std::endl;
		return false;
	}

	Smooth::AddProperties(mesh, mode);

    // main function
    Smooth::OriginalRank(mesh);
	Smooth::VertexSmooth(mesh, mode, numberOfIterations);

    // add vertex normals
	Smooth::VertexNormal(mesh);

	//Smooth.Save(argv[4]);
	std::cout << "Saving : " << outputFilename << '\n';
	if (!OpenMesh::IO::write_mesh(mesh, outputFilename))
	{
		std::cerr << "Error: cannot write mesh to " << outputFilename << std::endl;
		return false;
	}

    // print simple rank
	Smooth::SimpleRank(mesh, mode);

	//std::cout << "Press any key to continue";
	//std::cin.get();
    return 0;
}