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
  void SetBrickDimensions(unsigned int _xBrickDim,
                          unsigned int _yBrickDim,
                          unsigned int _zBrickDim);
  void SetPaddingWidth(unsigned int _paddingWidth);

  // Do everything!
  bool Construct();

private:
  Forge();
  Forge(const Forge&);

  std::string inFilename_;
  std::string outFilename_;
  BricksHeader * header_;
  std::vector<Brick<real>*> bricks_;

  unsigned int structure_;
  unsigned int xBrickDim_;
  unsigned int yBrickDim_;
  unsigned int zBrickDim_;
  unsigned int paddingWidth_;

  // TODO work in progress

  const std::string tempFilename_ = "octree.tmp";
  const std::string tspFilename_ = "/home/vsand/OpenSpace/output.tsp";


  // Read header from file and set brick header info
  bool CreateHeader();
  // Create an octree for every timestep, save in one common file
  bool CreateOctree();
  // Delete the created temp file
  bool DeleteTempFile();
  // Use temp octrees to create TSP tree
  bool ConstructTSPTree();

  std::fstream instream_;

  // Points to first data entry after header
  std::ios::pos_type headerOffset_;

  unsigned int nrBricksBaseLevel_;
  unsigned int nrLevels_;
  unsigned int nrBricksPerOctree_;

  // Calculate Z-order index from x, y, z coordinates
  uint32_t ZOrder(uint16_t x, uint16_t y, uint16_t z);



};


}

#endif

