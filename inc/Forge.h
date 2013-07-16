/*
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#ifndef FORGE_H_
#define FORGE_H_

#include <string>
#include <vector>
#include <fstream>

#define real float

namespace osp {
  
class BricksHeader;
template <class T>
class Brick;

class Forge {
public:
  static Forge * New();
  ~Forge();

  void SetInFilename(std::string _inFilename);
  void SetOutFilename(std::string _outFilename);
  void SetStructure(unsigned int _structure);
  void SetBrickDimensions(unsigned int _xBrickDim);
  void SetPaddingWidth(unsigned int _paddingWidth);

  // Do everything!
  bool Construct();

private:
  Forge();
  Forge(const Forge&);

  std::string inFilename_;
  std::string outFilename_;
  std::vector<Brick<real>*> bricks_;

  // Metadata to be read
  unsigned int structure_;
  unsigned int dataDimensionality_;
  unsigned int brickDim_;
  unsigned int numBricks_;
  unsigned int numTimesteps_;
  unsigned int dim_;
  unsigned int paddingWidth_;
  unsigned int dataSize_;

  // Additional metadata
  unsigned int nrBricksBaseLevel_;
  unsigned int nrLevels_;
  unsigned int nrBricksPerOctree_;
  unsigned int paddedDim_;
  unsigned int paddedBrickDim_;

  // TODO work in progress
  const std::string tempFilename_ = "octree.tmp";
  const std::string tspFilename_ = "/home/vsand/OpenSpace/output.tsp";

  // Read metadata
  bool ReadMetadata();
  // Create an octree for every timestep, save in one common file
  bool CreateOctree();
  // Delete the created temp file
  bool DeleteTempFile();
  // Use temp octrees to create TSP tree
  bool ConstructTSPTree();
  // Use temo octrees to create TSP tree (spatial ordering)
  bool ConstructTSPTreeSpatial();

  std::fstream instream_;

  // Points to first data entry after header
  std::ios::pos_type headerOffset_;
  // Calculate Z-order index from x, y, z coordinates
  uint32_t ZOrder(uint16_t x, uint16_t y, uint16_t z);


};


}

#endif

