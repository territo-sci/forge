/*
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 * Forge takes a VDF input file, rearranges and filters the data
 * to produce a corresponding TSP tree.
 *
 */

#ifndef FORGE_H_
#define FORGE_H_

// Make sure we get 64 bits for offset
#define _FILE_OFFSET_BITS 64
// For easy switching between offset types
#define off off64_t

#include <string>
#include <vector>
#include <list>

namespace osp {
  
class BricksHeader;
template <class T>
class Brick;

class Forge {
public:
  static Forge * New();
  ~Forge();

  // Setters for various variables
  void SetInFilename(std::string _inFilename);
  void SetOutFilename(std::string _outFilename);
  void SetBrickDimensions(unsigned int _xBrickDim,
                          unsigned int _yBrickDim,
                          unsigned int _zBrickDim);

  // Execute all construction steps
  bool Construct();

private:
  Forge();
  Forge(const Forge&);

  std::string inFilename_;
  std::string outFilename_;
  std::string tempFilename_;

  const unsigned int paddingWidth_ = 1;

  std::vector<Brick<float>*> bricks_;

  // Data that ends up in the out file
  unsigned int gridType_;
  unsigned int numTimesteps_;
  unsigned int xDim_;
  unsigned int yDim_;
  unsigned int zDim_;
  unsigned int xBrickDim_;
  unsigned int yBrickDim_;
  unsigned int zBrickDim_;
  unsigned int xNumBricks_;
  unsigned int yNumBricks_;
  unsigned int zNumBricks_;
  unsigned int dataSize_;

  // Additional metadata
  unsigned int numBricksBaseLevel_;
  unsigned int numLevels_;
  unsigned int numBricksPerOctree_;
  unsigned int numBricksTotal_;
  unsigned int xPaddedDim_;
  unsigned int yPaddedDim_;
  unsigned int zPaddedDim_;
  unsigned int xPaddedBrickDim_;
  unsigned int yPaddedBrickDim_;
  unsigned int zPaddedBrickDim_;

  // Read metadata from VDF file and calculate additional metadata
  bool ReadMetadata();
  // Create an octree for every timestep, save in one common file
  bool CreateOctree();
  // Delete the created temp files
  bool DeleteTempFiles();
  // Use temp octrees to create TSP tree 
  bool ConstructTSPTree();

  // Points to first data entry after header
  off headerOffset_;

  // Calculate Z-order index from x, y, z coordinates
  uint32_t ZOrder(uint16_t x, uint16_t y, uint16_t z);

};

}

#endif

