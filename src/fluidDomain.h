#pragma once

#include "Geometry.h"
#include "Mesh.h"

class FluidDomain
{
public:
	FluidDomain(Geometry* geometry, const int& index = 0);

	~FluidDomain();

	NodeMesh* getNode(const int& index);

	void addNode(const int& index, const bounded_vector<double,2>& coord);

	void addSurfaceMaterial(const std::vector<PlaneSurface*> surfaces, const double& young, const double& poisson, const double& density = 0.0, const std::string& behaviour = "PLANE_STRESS");

	void addElement(const int& index, const std::vector<int>& nodes, const int& materialIndex, const double& thickness, const std::string& elementType);

	void generateMesh(const std::string& elementType = "T3", const std::string& algorithm = "AUTO", std::string geofile = std::string(), const std::string& gmshPath = std::string(), const bool& plotMesh = true, const bool& showInfo = false);

	void readInput(const std::string& inputFile, const bool& deleteFiles = true);

public:
	int index_;
	Geometry * geometry_;
	std::vector<NodeMesh*> nodes_;
	std::vector<ElementMesh*> elements_;
	std::vector<Material*> materials_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

std::vector<std::string> split(std::string str, std::string delim)
{
	std::istringstream is(str);
	std::vector<std::string> values;
	std::string token;
	while (getline(is, token, ' '))
		values.push_back(token);
	return values;
}

std::vector<int> intersection(std::vector<int> vector1, std::vector<int> vector2)
{
	std::sort(vector1.begin(), vector1.end());
	std::sort(vector2.begin(), vector2.end());

	std::vector<int> vector3;
	std::set_intersection(vector1.begin(), vector1.end(), vector2.begin(), vector2.end(), std::back_inserter(vector3));

	return vector3;
}

std::vector<double> splitdouble(std::string str, std::string delim)
{
	std::istringstream is(str);
	std::vector<double> values;
	std::string token;
	while (getline(is, token, ' '))
		values.push_back(stod(token));
	return values;
}
std::vector<bool> splitbool(std::string str, std::string delim)
{
	std::istringstream is(str);
	std::vector<bool> values;
	std::string token;
	while (getline(is, token, ' '))
		values.push_back(token == "1");
	return values;
}

FluidDomain::FluidDomain(Geometry* geometry, const int& index)
{
	geometry_ = geometry;
	index_ = index;
}

FluidDomain::~FluidDomain() {}

NodeMesh* FluidDomain::getNode(const int& index)
{
	return nodes_[index];
}

void FluidDomain::addNode(const int& index, const bounded_vector<double,2>& coord)
{
	NodeMesh* n = new NodeMesh(index, coord);
	nodes_.push_back(n);
}

void FluidDomain::addElement(const int& index, const std::vector<int>& nodesIndex, const int& materialIndex, const double& thickness, const std::string& elementType)
{
	std::vector<NodeMesh*> nodes;
	nodes.reserve(nodesIndex.size());
	for (const int& i : nodesIndex)
		nodes.push_back(nodes_[i]);
	ElementMesh* e = new ElementMesh(index, nodes, materials_[materialIndex], thickness);
	elements_.push_back(e);
}

void FluidDomain::addSurfaceMaterial(const std::vector<PlaneSurface*> surfaces, const double& young, const double& poisson, const double& density, const std::string& behaviour)
{
	int index = materials_.size();
	Material* m = new Material(index, young, poisson, density, behaviour);
	materials_.push_back(m);
	for (PlaneSurface* surface : surfaces)
		surface->setMaterial(m);
}

void FluidDomain::readInput(const std::string& inputFile, const bool& deleteFiles)
{
	//defyning the maps that are used to store the elements information
	std::unordered_map<int, std::string> gmshElement = { {1, "line"}, {2, "triangle"}, {3, "quadrilateral"}, {8, "line3"}, {9, "triangle6"}, {10, "quadrilateral9"}, {15, "vertex"}, {16, "quadrilateral8"}, {20, "triangle9"}, {21, "triangle10"}, {26, "line4"}, {36, "quadrilateral16"}, {39, "quadrilateral12"} };
	std::unordered_map<std::string, int> numNodes = { {"vertex", 1}, {"line", 2}, {"triangle", 3}, {"quadrilateral", 4}, {"line3", 3}, {"triangle6", 6}, {"quadrilateral8", 8}, {"quadrilateral9", 9}, {"line4", 4}, {"triangle", 9}, {"triangle10", 10}, {"quadrilateral12", 12}, {"quadrilateral16", 16}};
	std::unordered_map<std::string, std::string> supportedElements = { {"triangle", "T3"}, {"triangle6", "T6"}, {"triangle10", "T10"}, {"quadrilateral", "Q4"}, {"quadrilateral8", "Q8"}, {"quadrilateral9", "Q9"}, {"quadrilateral12", "Q12"}, {"quadrilateral16", "Q16"} };
	std::unordered_map<Line*, std::vector< std::vector<int> >> lineElements;

	//opening the .msh file
	std::ifstream file(inputFile);
	std::string line;
	std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);
  

	//reading physical entities
	int number_physical_entities;
	file >> number_physical_entities;
	std::getline(file, line);
	std::unordered_map<int, std::string> physicalEntities;
	physicalEntities.reserve(number_physical_entities);

	std::cout << "NSDFSDAA 2 " << number_physical_entities << std::endl;

	for (int i = 0; i < number_physical_entities; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line, " ");
		int index;
		std::istringstream(tokens[1]) >> index;
		physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
	}
	std::getline(file, line); std::getline(file, line);
	//reading nodes
	int number_nodes;
	file >> number_nodes;
	nodes_.reserve(number_nodes);
	std::getline(file, line);
	for (int i = 0; i < number_nodes; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line, " ");
		bounded_vector<double,2> coord;
		std::istringstream(tokens[1]) >> coord(0);
		std::istringstream(tokens[2]) >> coord(1);
		addNode(i, coord);
	}
	std::getline(file, line); std::getline(file, line);
	//reading elements
	int number_elements;
	file >> number_elements;
	elements_.reserve(number_elements);
	std::getline(file, line);
	int cont = 0;
	for (int i = 0; i < number_elements; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line, " ");
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		std::string elementType = gmshElement[values[1]];
		int number_nodes_per_element = numNodes[elementType];
		std::vector<int> elementNodes;
		elementNodes.reserve(number_nodes_per_element);
		for (size_t j = 5 ; j < values.size(); j++)
			elementNodes.push_back(values[j]-1);
		std::string name = physicalEntities[values[3]];
		//Adding 2D elements to surfaces
		if (name[0] == 's')
		{
			if (supportedElements.find(elementType) == supportedElements.end())
			{
				std::cout << elementType << " is not supported.\n";
				exit(EXIT_FAILURE);
			}
			PlaneSurface* object = geometry_->getPlaneSurface(name);
			int materialIndex = object->getMaterial()->getIndex();
			double thickness = object->getThickness();
			addElement(cont, elementNodes, materialIndex, thickness, supportedElements[elementType]);
			object->addElementToSurface(elements_[cont]);
			//Verifying if an element side touches a line
			//for ( auto& pair : lineElements)
			//{
			//	for (std::vector<int> nodes : pair.second)
			//	{
			//		std::vector< std::vector<int> >::iterator it;
			//		it = std::find(pair.second.begin(), pair.second.end(), nodes);
			//		std::sort(nodes.begin(), nodes.end());
			//		if (intersection(nodes, elementNodes) == nodes)
			//		{
			//			pair.first->addElementsToLine(elements_[cont]);
			//			it = pair.second.erase(it);
			//			//if (pair.second.size() == 0)
			//				lineElements.erase(pair.first);
			//		}
			//	}
			//}

			for (auto it = lineElements.begin(); it != lineElements.end();)
			{
				auto& key = it->first;
				auto& value = it->second;
				int it2 = 0;
				for (auto it1 = value.begin(); it1 != value.end();)
				{
					auto nodes = lineElements[key][it2];
					std::sort(nodes.begin(), nodes.end());
					if (intersection(nodes, elementNodes) == nodes)
					{
						key->addElementsToLine(elements_[cont]);
						it1 = value.erase(it1);
					}
					else
					{
						it1++;
						it2++;
					}
				}
				if (value.size() == 0)
					it = lineElements.erase(it);
				else
					it++;

			}
			
			cont += 1;
		}
		//Adding 1D elements to lines
		else if (name[0] == 'l')
		{
			Line* object = geometry_->getLine(name);
			std::vector<NodeMesh*> nodes;
			nodes.reserve(elementNodes.size());
			for (int index : elementNodes)
				nodes.push_back(nodes_[index]);
			object->addNodesToLine(nodes);
			if (lineElements.count(object) != 1)
				lineElements[object] = std::vector< std::vector<int> >();
			lineElements[object].push_back(elementNodes);
		}
		//Adding a node to point
		else
		{
			Point* object = geometry_->getPoint(name);
			object->addNodeToPoint(getNode(elementNodes[0]));
		}
	}
	//Closing the file
	file.close();
	if (deleteFiles)
		system((remove2 + inputFile).c_str());
	return;
}

void FluidDomain::generateMesh(const std::string& elementType, const std::string& algorithm, std::string geofile, const std::string& gmshPath, const bool& plotMesh, const bool& showInfo)
{
	std::pair<std::string, bool> pair = createMesh(geometry_, elementType, algorithm, geofile, gmshPath, plotMesh, showInfo);

	//readInput(pair.first, pair.second);
}