#pragma once

#include <fstream>
#include <stdio.h>
#include "Geometry.h"

#ifdef _WIN32
#include <direct.h>
#define getCurrentDir _getcwd
//#define remove "del "
#else
#include <unistd.h>
#define getCurrentDir getcwd
//#define remove "rm "
#endif

std::string remove2 = "rm ";
//std::string getCurrentDir = "getcwd";
std::string getCurrentWorkingDir();

std::pair<std::string, bool> createMesh(Geometry* geometry, const std::string& elementType = "T3", const std::string& algorithm = "AUTO", std::string geofile = std::string(), const std::string& gmshPath = std::string(), const bool& plotMesh = true, const bool& showInfo = false);


///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------


std::string getCurrentWorkingDir()
{
	char buff[FILENAME_MAX];
	getCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

std::pair<std::string, bool> createMesh(Geometry* geometry, const std::string& elementType, const std::string& algorithm, std::string geofile, const std::string& gmshPath, const bool& plotMesh, const bool& showInfo)
{
	std::string mshfile;
	if (elementType[0] == 'Q' || elementType[0] == 'q')
	{
		std::stringstream text;
		text << "Recombine Surface {";
		/*for (int i = 0; i < geometry->getNumberOfPlaneSurfaces(); i++)
		{
			text << geometry->getPlaneSurface(i)->getName();
			if (i != (geometry->getNumberOfPlaneSurfaces() - 1))
				text << ", ";
		}*/
		int cont = 0;
		for (const auto& pair : geometry->getPlaneSurfaces())
		{
			text << pair.first;
			if (cont != (geometry->getNumberOfPlaneSurfaces() - 1))
				text << ", ";
			cont += 1;
		}
		text << "};";
		geometry->appendGmshCode(text.str());
	}
	bool deleteFiles = geofile.empty();
	if (deleteFiles)
	{
		mshfile = geofile + "temp.msh";
		geofile = geofile + "temp.geo";
	}
	else {
		mshfile = geofile + ".msh";
		geofile = geofile + ".geo";
	}
	std::ofstream file(geofile);
	file << geometry->getGmshCode();
	file.close();

	std::string gmshExe = (gmshPath.empty()) ? getCurrentWorkingDir() + "/src/gmsh" : gmshPath;
	std::string cmd = gmshExe;
	cmd += " -2 -clscale 1.0 " + geofile + " -o " + mshfile;

	if (algorithm == "FRONT")
		cmd += " -algo front2d";
	else if (algorithm == "DELAUNAY")
		cmd += " -algo del2d";
	else if (algorithm == "ADAPT")
		cmd += " -algo meshadapt";
	else if (algorithm == "PACK")
		cmd += " -algo pack";
	else if (algorithm == "QUAD")
		cmd += " -algo delquad";
	else if (algorithm == "BAMG")
		cmd += " -algo bamg";

	if (elementType == "T3" || elementType == "Q4")
		cmd += " -order 1";
	else if (elementType == "T6" || elementType == "Q8" || elementType == "Q9") {
		cmd += " -order 2";
		if (elementType == "Q8")
			cmd += " -string Mesh.SecondOrderIncomplete=1;";
	}
	else if (elementType == "T10" || elementType == "Q12" || elementType == "Q16") {
		cmd += " -order 3";
		if (elementType == "Q12")
			cmd += " -string Mesh.SecondOrderIncomplete=1;";
	}
	else
	{
		std::cout << "Element " << elementType << " is not supported.";
		exit(EXIT_FAILURE);
	}

	if (!showInfo) cmd += " -v 0";

	system(cmd.c_str());

	if (deleteFiles)
		system((remove2 + geofile).c_str());

	if (plotMesh)
		system((gmshExe + " " + mshfile).c_str());

	return std::make_pair(mshfile, deleteFiles);
}