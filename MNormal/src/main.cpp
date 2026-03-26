/**
 * @file main.cpp
 * @author yhf (hfye@zju.edu.cn)
 * @brief this is a mesh tool modify surface mesh by virtual topologic generation so that enable to generate mutiple normals in boundary layer mesh generation.  
 * @version 1.6
 * @date 2022-01-24
 * 
 * @copyright Copyright (c) 2022
 * 
 */

 //  $$\      $$\                     $$\        $$$$$$\  $$\                                $$$$$$\                     
 //  $$$\    $$$ |                    $$ |      $$  __$$\ $$ |                              $$  __$$\                    
 //  $$$$\  $$$$ | $$$$$$\   $$$$$$$\ $$$$$$$\  $$ /  \__|$$$$$$$\   $$$$$$\  $$$$$$\$$$$\  $$ /  \__|$$$$$$\   $$$$$$\  
 //  $$\$$\$$ $$ |$$  __$$\ $$  _____|$$  __$$\ $$ |      $$  __$$\  \____$$\ $$  _$$  _$$\ $$$$\    $$  __$$\ $$  __$$\ 
 //  $$ \$$$  $$ |$$$$$$$$ |\$$$$$$\  $$ |  $$ |$$ |      $$ |  $$ | $$$$$$$ |$$ / $$ / $$ |$$  _|   $$$$$$$$ |$$ |  \__|
 //  $$ |\$  /$$ |$$   ____| \____$$\ $$ |  $$ |$$ |  $$\ $$ |  $$ |$$  __$$ |$$ | $$ | $$ |$$ |     $$   ____|$$ |      
 //  $$ | \_/ $$ |\$$$$$$$\ $$$$$$$  |$$ |  $$ |\$$$$$$  |$$ |  $$ |\$$$$$$$ |$$ | $$ | $$ |$$ |     \$$$$$$$\ $$ |      
 //  \__|     \__| \_______|\_______/ \__|  \__| \______/ \__|  \__| \_______|\__| \__| \__|\__|      \_______|\__|      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
 //                                                                                                                      
#include<iostream>
#include<sstream>
#include "MNormalMesh.h"
#include "../include/common.h"
#include "ChamferBehavior.h"
#include "../third/cli11/CLI11.hpp"
#include <spdlog/spdlog.h>
#define F16
int main(int argc, char** argv) {


		
	// read logo from file

		
	ChamferBehavior behavior;
	CLI::App app{ "MeshChamfer"};
	app.description("A plane mesh is processed to make it suitable for multi-normal boundary mesh generation.");
	app.add_option("-i", behavior.input_filename, "input filename. (string, required, supported format: vtk, mesh, pls, obj)")->required();
	


	
	/// prase the input command by CLI
	try {
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError &e) {
		return app.exit(e);
	}
	cout<<"================================================================"<<endl;
	cout<<GREEN<<"Mesh Chamfer, a solution for mutiple normal generation."<<RESET<<std::endl;
	cout <<GREEN<< "Input Name:  	"<<RESET << behavior.input_filename << endl;
	cout <<GREEN<< "Output Name:	" <<RESET<< behavior.getOutputName(".pls") << endl;
	cout<<"================================================================"<<endl;


	spdlog::info("Reading the model!");
	
	/// read input and generate initial visibility graph
	MNormalMesh chamfer; // create one chamfer 
	chamfer.SetBehavior(behavior);
	chamfer.ReadPls();
	spdlog::info("Done!");

	chamfer.CalculateMultiNormal();
	int faceCount;
	chamfer.BuildTopo(faceCount);
	
	spdlog::info("Handling output mesh!");
	chamfer.WritePls();
	chamfer.WriteVtk();
	chamfer.WriteNorm();
	chamfer.GenerateFirstLayer(0.01);

	spdlog::info("Job Finished.");
}