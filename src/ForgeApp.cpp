/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#include <Forge.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace osp;

int main() {


  // Very simple and short config file
  // TODO nicer implementation
  std::string config = "config/forgeConfig.txt";
  std::ifstream in;
  in.open(config.c_str(), std::ifstream::in);
  if (!in.is_open()) {
    std::cout << "Could not open " << config << std::endl;
    exit(1);
  }

  std::string line;
  std::string inFilename = "notSet";
  std::string outFilename = "notSet";
  unsigned int xBrickDim = 0;
  unsigned int yBrickDim = 0;
  unsigned int zBrickDim = 0;
  float spatialScaling = 0.f;
  float temporalScaling = 0.f;
  while (std::getline(in, line)) {
    // Ignore empty lines and comments
    if (!line.empty() && line.at(0) != '#') {
      // Read variable name
      std::stringstream ss;
      ss.str(line);
      std::string var;
      ss >> var;
      // Read value
      if (var == "in_filename") {
        ss >> inFilename;
      } else if (var == "out_filename") {
        ss >> outFilename;
      } else if (var == "brick_dimensions") {
        ss >> xBrickDim;
        ss >> yBrickDim;
        ss >> zBrickDim;
      } else if (var == "spatial_scaling") {
        ss >> spatialScaling;
      } else if (var == "temporal_scaling") {
        ss >> temporalScaling;
      } else {
        std::cout << "Variable " << var << " not recognized" << std::endl;
        exit(1);
      }
    }
  }

  in.close();

  std::cout << "In filename: " << inFilename << std::endl;
  std::cout << "Out filename: " << outFilename << std::endl;
  std::cout << "Brick dimensions: " << xBrickDim << " " << yBrickDim 
    << " " << zBrickDim << std::endl;
  std::cout << "Spatial scaling: " << spatialScaling << std::endl;
  std::cout << "Temporal scaling: " << temporalScaling << std::endl;

  Forge* forge = Forge::New();
  
  forge->SetInFilename(inFilename);
  forge->SetOutFilename(outFilename);
  forge->SetBrickDimensions(xBrickDim, yBrickDim, zBrickDim);
  forge->SetSpatialScaling(spatialScaling);
  forge->SetTemporalScaling(temporalScaling);
 
  // Construct TSP tree and calculate errors 
  if (!forge->Construct()) { 
    std::cerr << "Failed to construct TSP tree" << std::endl;
    exit(1);
  }

  delete forge;
  exit(0);
}
