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

#include <string>
#include <vector>
#include <list>

#include <glm/glm.hpp>

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
  void SetBrickDimensions(size_t _xBrickDim,
                          size_t _yBrickDim,
                          size_t _zBrickDim);

  // Execute all construction steps
  bool Construct();

private:
  Forge();
  Forge(const Forge&);

  std::string inFilename_;
  std::string outFilename_;
  std::string tempFilename_;

  const size_t paddingWidth_ = 1;

  std::vector<Brick<float>*> bricks_;

  // Data that ends up in the out file
  size_t gridType_;
  size_t numTimesteps_;
  size_t xDim_;
  size_t yDim_;
  size_t zDim_;
  size_t xBrickDim_;
  size_t yBrickDim_;
  size_t zBrickDim_;
  size_t xNumBricks_;
  size_t yNumBricks_;
  size_t zNumBricks_;
  size_t dataSize_;

  // Additional metadata
  size_t numBricksBaseLevel_;
  size_t numLevels_;
  size_t numBricksPerOctree_;
  size_t numBricksTotal_;
  size_t xPaddedDim_;
  size_t yPaddedDim_;
  size_t zPaddedDim_;
  size_t xPaddedBrickDim_;
  size_t yPaddedBrickDim_;
  size_t zPaddedBrickDim_;

  // Read metadata from VDF file and calculate additional metadata
  bool ReadMetadata();
  // Create an octree for every timestep, save in one common file
  bool CreateOctrees();
  // Delete the created temp files
  bool DeleteTempFiles();
  // Use temp octrees to create TSP tree 
  bool ConstructTSPTree();

  // Points to first data entry after header
  size_t headerOffset_;

  bool buildDataLevels(std::FILE* file, size_t timestep, std::vector< std::vector<float> >& levelData);
  bool createPadding(std::vector< std::vector<float> >& levelData, std::vector< std::vector<float> >& paddedLevelData);
  bool buildOctree(std::vector< std::vector<float> >& paddedLevelData, std::vector< Brick<float>* >& octreeBricks);

  // Calculate Z-order index from x, y, z coordinates
  uint32_t ZOrder(uint16_t x, uint16_t y, uint16_t z);

  glm::ivec3 linearToCartesian(size_t linearCoords, int dim);
  glm::ivec3 linearToCartesian(size_t linearCoords, int xDim, int yDim, int zDim);
  glm::ivec3 linearToCartesian(size_t linearCoords, glm::ivec3 dim);
  size_t cartesianToLinear(glm::ivec3 cartesianCoords, int dim);
  size_t cartesianToLinear(glm::ivec3 cartesianCoords, int xDim, int yDim, int zDim);
  size_t cartesianToLinear(glm::ivec3 cartesianCoords, glm::ivec3 dim);

};

} // namespace osp

#endif

