/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 * Note: Uses C-style I/O to handle large files. This functionality
 * might be implementation/platform dependent!
 *
 * Using fseeko/ftello instead of fseek/ftell to ensure 64 bit offsets.
 *
 */

#include <Forge.h>
#include <Brick.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <array>
#include <boost/filesystem.hpp>
#include <stdio.h>

namespace osp {

Forge * Forge::New() {
  return new Forge();
}

Forge::Forge() 
  : inFilename_("NotSet"), 
    outFilename_("NotSet"), 
    tempFilename_("temp.tmp"),
    gridType_(0), 
    numTimesteps_(0),
    xDim_(0),
    yDim_(0),
    zDim_(0),
    xBrickDim_(0),
    yBrickDim_(0),
    zBrickDim_(0),
    xNumBricks_(0),
    yNumBricks_(0),
    zNumBricks_(0),
    dataSize_(0),
    numBricksBaseLevel_(0),
    numLevels_(0),
    numBricksPerOctree_(0),
    numBricksTotal_(0),
    xPaddedDim_(0),
    yPaddedDim_(0),
    zPaddedDim_(0),
    xPaddedBrickDim_(0),
    yPaddedBrickDim_(0),
    zPaddedBrickDim_(0) {
}

Forge::~Forge() {
  for (auto it=bricks_.begin(); it!=bricks_.end(); ++it) {
    delete *it;
  }
}

void Forge::SetInFilename(std::string _inFilename) {
  inFilename_ = _inFilename;
}

void Forge::SetOutFilename(std::string _outFilename) {
  outFilename_ = _outFilename;
}

void Forge::SetBrickDimensions(unsigned int _xBrickDim,
                               unsigned int _yBrickDim,
                               unsigned int _zBrickDim) {
  xBrickDim_ = _xBrickDim;
  yBrickDim_ = _yBrickDim;
  zBrickDim_ = _zBrickDim;
}

bool Forge::Construct() {
  
  // Read metadata from VDF file
  if (!ReadMetadata()) {
    std::cerr << "Error: Could not read metadata" << std::endl;
    return false;
  }

  // Create one octree per timestep, save in temporary file
  if (!CreateOctrees()) {
    std::cerr << "Error: Failed to create temp octree file" << std::endl;
    return false;
  }

  // Construct the TSP tree from the temporary files
  if (!ConstructTSPTree()) {
    std::cerr << "Error: Failed to construct TSP tree" << std::endl;
    return false;
  }

  // Delete temporary files
  if (!DeleteTempFiles()) {
    std::cerr << "Failed to delete temp files, but that's okay" << std::endl;
    std::cerr << "Make sure to clean up the Forge directory manually" <<
      std::endl;
  }

  return true;
}


bool Forge::ReadMetadata() {

  std::cout << "Reading metadata" << std::endl;

  if (xBrickDim_ == 0 || yBrickDim_ == 0 || zBrickDim_ == 0) {
    std::cerr << "One or more brick dimension is zero!"<< std::endl;
    return false;
  }

  // TODO this is pretty pointless, can be removed
  dataSize_ = sizeof(float);

  // Read  header data from VDF file
  std::FILE *in = fopen(inFilename_.c_str(), "r");
  if (!in) {
    std::cerr << "Failed to open " << inFilename_ << std::endl;
    return false;
  }
  
  size_t s = sizeof(unsigned int);
  fread(reinterpret_cast<void*>(&gridType_), s, 1, in);
  fread(reinterpret_cast<void*>(&numTimesteps_), s, 1, in);
  fread(reinterpret_cast<void*>(&xDim_), s, 1, in);
  fread(reinterpret_cast<void*>(&yDim_), s, 1, in);
  fread(reinterpret_cast<void*>(&zDim_), s, 1, in);

  // TODO support non-full BST trees. Right now, the number of timesteps
  // needs to be a power of two. Abort if it's not.

  if (xDim_ % xBrickDim_ != 0 ||
      yDim_ % yBrickDim_ != 0 ||
      zDim_ % zBrickDim_ != 0) {
    std::cerr << "Error: Voxel/brick dimension mismatch!" << std::endl;
    fclose(in);
    return false;
  }

  // Number of bricks per axis
  xNumBricks_ = xDim_ / xBrickDim_;
  yNumBricks_ = yDim_ / yBrickDim_;
  zNumBricks_ = zDim_ / zBrickDim_;
  // Dimensions of padded volume 
  xPaddedDim_ = xDim_ + paddingWidth_*2;
  yPaddedDim_ = yDim_ + paddingWidth_*2;
  zPaddedDim_ = zDim_ + paddingWidth_*2;
  // Dimensions of padded bricks
  xPaddedBrickDim_ = xBrickDim_ + paddingWidth_*2;
  yPaddedBrickDim_ = yBrickDim_ + paddingWidth_*2;
  zPaddedBrickDim_ = zBrickDim_ + paddingWidth_*2;

  std::cout << std::endl << "FORGE METADATA" << std::endl;
  std::cout << "Read from " << inFilename_ << " complete!" << std::endl;
  std::cout << "Grid type: " << gridType_ << std::endl;
  std::cout << "Number of timesteps: " << numTimesteps_ << std::endl;
  std::cout << "Dimensions: " << xDim_ << " x " << yDim_ <<  
               " x " << zDim_ << std::endl;
  std::cout << "Brick dimensions: " << xBrickDim_ << " x " << yBrickDim_ <<
               " x " << zBrickDim_ << std::endl;
  std::cout << "Padded dimensions: " << xPaddedDim_ << " x " << yPaddedDim_ <<
               " x " << zPaddedDim_ << std::endl;
  std::cout << "Padded brick dimensions: " << xPaddedBrickDim_ << " x " <<
                yPaddedBrickDim_ << " x " << zPaddedBrickDim_ << std::endl;
  std::cout << "Number of bricks: " << xNumBricks_ << " x " << yNumBricks_ << 
               " x " << zNumBricks_ << std::endl;
  std::cout << "Data size (bytes): " << dataSize_ << std::endl;
  std::cout << "Out file name: " << outFilename_ << std::endl;

  // Calculate some common things

  // Number of bricks in the base (leaf) level
  numBricksBaseLevel_ = xNumBricks_*yNumBricks_*zNumBricks_;
  // Number of octree levels
  // TODO support different axis dimensions
  numLevels_ = log(xNumBricks_)/log(2) + 1;
  // Number of bricks per octree (used for offsets)
  numBricksPerOctree_ = (pow(8, numLevels_) - 1) / 7;
  // Number of bricks total
  numBricksTotal_ = numBricksPerOctree_ * (2*numTimesteps_ - 1);

  std::cout << "Number of bricks in base octree level: " 
    << numBricksBaseLevel_ << std::endl;
  std::cout << "Number of levels in octree: " << numLevels_ << std::endl;
  std::cout << "Number of bricks in octree: " 
    << numBricksPerOctree_ << std::endl;
  std::cout << "Number of bricks total " << numBricksTotal_ << std::endl;

  // Save position of first data entry after header
  headerOffset_ = ftello(in);
  
  fclose(in);

  return true;
}

bool Forge::CreateOctrees() {

  // Init in and out files

  std::FILE *out = fopen(tempFilename_.c_str(), "w");
  if (!out) {
    std::cerr << "Failed to init " << tempFilename_ << std::endl;
    return false;
  }

  std::FILE *in = fopen(inFilename_.c_str(), "r");
  if (!in) {
    std::cerr << "Failed to init " << inFilename_ << std::endl;
    fclose(out);
    return false;
  }

  // Loop over all timesteps 
  for (unsigned int timestep = 0; timestep < numTimesteps_; timestep++) {
    std::cout << "Constructing octree for timestep " << timestep << "/" << numTimesteps_ << "\r" << std::flush;

    unsigned int nOctreeLevels = log(numBricksBaseLevel_)/log(8) + 1;
    std::vector< std::vector<float> > levelData(nOctreeLevels);

    if (!buildDataLevels(in, timestep, levelData)) {
      std::cerr << std::endl << "Failed to build data levels" << std::endl;
      return false;
    }

    std::vector< std::vector<float> > paddedLevelData(nOctreeLevels);
    if (!createPadding(levelData, paddedLevelData)) {
      std::cerr << std::endl << "Failed to create padding" << std::endl;
      return false;
    }

    std::vector< Brick<float>* > octreeBricks(numBricksPerOctree_);
    if (!buildOctree(paddedLevelData, octreeBricks)) {
      std::cerr << std::endl << "Failed to build octree" << std::endl;
    }

    // Write octree to file
    for (auto it = octreeBricks.begin(); it != octreeBricks.end(); it++) {
      fwrite(reinterpret_cast<void*>(&(*it)->data_[0]),
             static_cast<size_t>((*it)->Size()), 1, out);
      delete *it;
      *it = nullptr;
    }
  }

  std::cout << "Constructing octree for timestep " << numTimesteps_ << "/" << numTimesteps_ << std::endl;

  fclose(in);
  fclose(out);

  return true;
}


bool Forge::buildDataLevels(std::FILE* file, unsigned int timestep, std::vector< std::vector<float> >& levelData) {
  unsigned int nOctreeLevels = levelData.size();

  unsigned int nBaseLevelVoxels = xDim_*yDim_*zDim_;
  levelData[0].resize(nBaseLevelVoxels);

  size_t timestepSize = nBaseLevelVoxels * dataSize_;
  size_t timestepOffset = static_cast<size_t>(timestep) * timestepSize + headerOffset_;
  fseeko(file, timestepOffset, SEEK_SET);
  fread(reinterpret_cast<void*>(&levelData[0][0]), timestepSize, 1, file);

  glm::ivec3 levelDim(xDim_, yDim_, zDim_);
  for (unsigned int level = 1; level < nOctreeLevels; level++) {
    unsigned int childLevel = level - 1;
    glm::ivec3 childLevelDim = levelDim;
    levelDim /= 2;
    unsigned int nLevelVoxels = levelDim.x * levelDim.y * levelDim.z;
    levelData[level].resize(nLevelVoxels);

    for (int z = 0; z < levelDim.z; z++) {
      for (int y = 0; y < levelDim.y; y++) {
        for (int x = 0; x < levelDim.x; x++) {
          glm::ivec3 voxelPos(x, y, z);

          float voxelData = 0.0;
          for (int childZ = 0; childZ < 2; childZ++) {
            for (int childY = 0; childY < 2; childY++) {
              for (int childX = 0; childX < 2; childX++) {
                glm::ivec3 localChildPos(childX, childY, childZ);
                glm::ivec3 childPos = voxelPos * 2 + localChildPos;
                unsigned int childVoxelIndex = cartesianToLinear(childPos, childLevelDim);
                voxelData += levelData[childLevel][childVoxelIndex];
              }
            }
          }
          voxelData /= 8.0;

          unsigned int voxelIndex = cartesianToLinear(voxelPos, levelDim);
          levelData[level][voxelIndex] = voxelData;
        }
      }
    }
  }

  return true;
}

bool Forge::createPadding(std::vector< std::vector<float> >& levelData, std::vector< std::vector<float> >& paddedLevelData) {
  unsigned int nOctreeLevels = log(numBricksBaseLevel_)/log(8) + 1;

  if (levelData.size() != nOctreeLevels || paddedLevelData.size() != nOctreeLevels) {
    return false;
  }

  glm::ivec3 paddingOffset(paddingWidth_);
  glm::ivec3 levelDim(xDim_, yDim_, zDim_);

  for (unsigned int level = 0; level < nOctreeLevels; level++) {
    glm::ivec3 paddedLevelDim = levelDim + paddingOffset * 2;
    unsigned int nPaddedLevelVoxels = paddedLevelDim.x * paddedLevelDim.y * paddedLevelDim.z;
    paddedLevelData[level].resize(nPaddedLevelVoxels);

    for (int z = 0; z < paddedLevelDim.z; z++) {
      for (int y = 0; y < paddedLevelDim.y; y++) {
        for (int x = 0; x < paddedLevelDim.x; x++) {
          glm::ivec3 paddedLevelVoxelPos(x, y, z);
          unsigned int paddedLevelVoxel = cartesianToLinear(paddedLevelVoxelPos, paddedLevelDim);

          glm::ivec3 samplePos = paddedLevelVoxelPos - paddingOffset;
          samplePos.x = glm::clamp(samplePos.x, 0, levelDim.x - 1);
          samplePos.y = glm::clamp(samplePos.y, 0, levelDim.y - 1);
          samplePos.z = glm::clamp(samplePos.z, 0, levelDim.z - 1);

          unsigned int levelVoxel = cartesianToLinear(samplePos, levelDim);
          paddedLevelData[level][paddedLevelVoxel] = levelData[level][levelVoxel];
        }
      }
    }

    levelDim /= 2;
  }

  return true;
}

bool Forge::buildOctree(std::vector< std::vector<float> >& paddedLevelData, std::vector< Brick<float>* >& octreeBricks) {
  unsigned int maxLevel = paddedLevelData.size() - 1;

  glm::ivec3 paddingOffset(paddingWidth_);
  glm::ivec3 brickDim(xBrickDim_, yBrickDim_, zBrickDim_);
  glm::ivec3 paddedBrickDim = brickDim + 2 * paddingOffset;

  unsigned int levelOffset = 0;

  for (unsigned int level = 0; level <= maxLevel; level++) {
    unsigned int depth = maxLevel - level;
    unsigned int nBricksPerDim = pow(2, depth);

    glm::ivec3 paddedLevelDataDim = brickDim;
    paddedLevelDataDim.x *= nBricksPerDim;
    paddedLevelDataDim.y *= nBricksPerDim;
    paddedLevelDataDim.z *= nBricksPerDim;
    paddedLevelDataDim += paddingOffset * 2;

    for (int brickZ = 0; brickZ < nBricksPerDim; brickZ++) {
      for (int brickY = 0; brickY < nBricksPerDim; brickY++) {
        for (int brickX = 0; brickX < nBricksPerDim; brickX++) {
          glm::ivec3 brickPos(brickX, brickY, brickZ);
          glm::ivec3 brickOffset = brickPos * brickDim;

          Brick<float> *brick = Brick<float>::New(paddedBrickDim.x, paddedBrickDim.y, paddedBrickDim.z, static_cast<float>(0));

          for (int voxelZ = 0; voxelZ < paddedBrickDim.z; voxelZ++) {
            for (int voxelY = 0; voxelY < paddedBrickDim.y; voxelY++) {
              for (int voxelX = 0; voxelX < paddedBrickDim.x; voxelX++) {
                glm::ivec3 voxelPos(voxelX, voxelY, voxelZ);
                glm::ivec3 voxelOffset = brickOffset + voxelPos;
                unsigned int voxel = cartesianToLinear(voxelOffset, paddedLevelDataDim);
                brick->SetData(voxelPos.x, voxelPos.y, voxelPos.z, paddedLevelData[level][voxel]);
              }
            }
          }

          unsigned int zOrderIdx = static_cast<unsigned int>(ZOrder(brickX, brickY, brickZ));
          unsigned int brickIndex = levelOffset + zOrderIdx;
          octreeBricks[brickIndex] = brick;
        }
      }
    }

    unsigned int nBricksInLevel = nBricksPerDim*nBricksPerDim*nBricksPerDim;
    levelOffset += nBricksInLevel;
  }

  return true;
}

glm::ivec3 Forge::linearToCartesian(unsigned int linearCoords, int dim) {
  return linearToCartesian(linearCoords, dim, dim, dim);
}

glm::ivec3 Forge::linearToCartesian(unsigned int linearCoords, int xDim, int yDim, int zDim) {
  return linearToCartesian(linearCoords, glm::ivec3(xDim, yDim, zDim));
}

glm::ivec3 Forge::linearToCartesian(unsigned int linearCoords, glm::ivec3 dim) {
  int z = linearCoords / (dim.x*dim.y);
  int y = (linearCoords / dim.x) % dim.y;
  int x = linearCoords % dim.x;
  return glm::ivec3(x, y, z);
}

unsigned int Forge::cartesianToLinear(glm::ivec3 cartesianCoords, int dim) {
  return cartesianToLinear(cartesianCoords, dim, dim, dim);
}

unsigned int Forge::cartesianToLinear(glm::ivec3 cartesianCoords, int xDim, int yDim, int zDim) {
  return cartesianToLinear(cartesianCoords, glm::ivec3(xDim, yDim, zDim));
}

unsigned int Forge::cartesianToLinear(glm::ivec3 cartesianCoords, glm::ivec3 dim) {
  return cartesianCoords.x + cartesianCoords.y*dim.x + cartesianCoords.z*dim.x*dim.y;
}

// Adapted from  http://graphics.stanford.edu/~seander/bithacks.htm
uint32_t Forge::ZOrder(uint16_t xPos, uint16_t yPos, uint16_t zPos) {
  uint32_t x = static_cast<uint32_t>(xPos);
  uint32_t y = static_cast<uint32_t>(yPos);
  uint32_t z = static_cast<uint32_t>(zPos);
  x = (x | (x << 16)) & 0x030000FF;
  x = (x | (x <<  8)) & 0x0300F00F;
  x = (x | (x <<  4)) & 0x030C30C3;
  x = (x | (x <<  2)) & 0x09249249;
  y = (y | (y << 16)) & 0x030000FF;
  y = (y | (y <<  8)) & 0x0300F00F;
  y = (y | (y <<  4)) & 0x030C30C3;
  y = (y | (y <<  2)) & 0x09249249;
  z = (z | (z << 16)) & 0x030000FF;
  z = (z | (z <<  8)) & 0x0300F00F;
  z = (z | (z <<  4)) & 0x030C30C3;
  z = (z | (z <<  2)) & 0x09249249;
  const uint32_t result = x | (y << 1) | (z << 2);
  return result;
}


bool Forge::DeleteTempFiles() {

  if (boost::filesystem::exists(tempFilename_)) {
    boost::filesystem::remove_all(tempFilename_);
  } else {
    std::cout <<"Warning: " << tempFilename_ << " does not exist" << std::endl;
  }

  unsigned int numBSTLevels = log(numTimesteps_)/log(2) + 1;
  for (unsigned int level=0; level<numBSTLevels; ++level) {

    std::stringstream ss;
    ss << level;
    std::string file = tempFilename_ + "." + ss.str() + ".tmp";

    if (boost::filesystem::exists(file)) {
      boost::filesystem::remove_all(file);
    } else {
      std::cout <<"Warning: " << file << " does not exist" << std::endl;
    }
  }

  return true;
}

bool Forge::ConstructTSPTree() {

  // Make sure temporary file exists
  if (!boost::filesystem::exists(tempFilename_)) {
    std::cout << "Error: temp file "<<tempFilename_<<" missing" << std::endl;
    return false;
  }

  // Numbers to keep track of
  unsigned int numOTNodes = numBricksPerOctree_;
  unsigned int numBrickVals = 
    xPaddedBrickDim_*yPaddedBrickDim_*zPaddedBrickDim_;

  std::cout << "Num nodes per OT: " << numOTNodes << std::endl;
  std::cout << "Num values per brick: " << numBrickVals << std::endl;

  // Keep track of the original number of timesteps in case we need to
  // adjust it. Both the original and the adjusted versions go in the header.
  unsigned int origNumTimesteps = numTimesteps_;

  // If the number of timesteps is not a power of two, we copy the last 
  // timestep to make the base level contain a number that is the next
  // power of two higher than the original number of timesteps
  bool powTwo = (numTimesteps_ & (numTimesteps_-1)) == 0;
  if (!powTwo) {

    // Find the next power of two higher than the number of timesteps
    unsigned int newNumTimesteps =  static_cast<unsigned int>(
      powf(2.f, ceilf(logf(static_cast<float>(numTimesteps_))/logf(2.f))));
    unsigned int timesToCopy = newNumTimesteps - numTimesteps_;

    std::cout << "Number of timesteps is not a power of two " << std::endl;
    std::cout << "Copying the last timestep " << timesToCopy << " times" <<
      " to get to " << newNumTimesteps << " timesteps " << std::endl; 

    // Read the last timestep
    size_t timestepSize = numBricksPerOctree_ * numBrickVals * sizeof(float);
    std::cout << "timestepSize: " << timestepSize << std::endl;
    size_t offset = static_cast<size_t>((numTimesteps_-1)*timestepSize);
    std::cout << "offset: " << offset << std::endl;

    std::FILE *in = fopen(tempFilename_.c_str(), "r");

    fseeko(in, 0, SEEK_END);
    size_t fileSize = ftello(in);
    std::cout << "File size: " << fileSize << std::endl;
    std::cout << "File size - timestep size: " << fileSize-timestepSize << std::endl;
    fseeko(in, offset, SEEK_SET);

    std::vector<float> buffer(timestepSize/sizeof(float));

    fread(reinterpret_cast<void*>(&buffer[0]), timestepSize, 1, in);
    fclose(in);

    std::FILE *out = fopen(tempFilename_.c_str(), "a+");
    if (!out) {
      std::cerr << "Failed to init " << tempFilename_ << std::endl;
      return false;
    }

    // Write the extra timesteps at the end of the file
    fseeko(out, 0, SEEK_END);
    for (unsigned int i=0; i<timesToCopy; ++i) {
      fwrite(reinterpret_cast<void*>(&buffer[0]), timestepSize, 1, out);  
    }

    fclose(out);

    // Use new number in calculations
    numTimesteps_ = newNumTimesteps;

  } else {
    std::cout << "Num timesteps is a power of two. Nice!" << std::endl;
  }

  // More numbers to keep track of
  unsigned int numBSTNodes = 2*numTimesteps_ - 1;
  unsigned int numBSTLevels = log(numTimesteps_)/log(2) + 1;

  std::cout << "Num nodes per BST: " << numBSTNodes << std::endl;
  std::cout << "Num BST levels: " << numBSTLevels << std::endl;

  // Append the temp file to match leaf BST level
  unsigned int BSTLevel = numBSTLevels - 1;

  std::string newFilename;

  { // Scoping stringstream
    std::stringstream ss;
    ss << BSTLevel;
    newFilename = tempFilename_ + "." + ss.str() + ".tmp";
  }

  std::cout << "Creating base level " << std::endl;

  // Create base level temp file by reversing the level order!
  
  { // Scoping files 

    std::FILE *in = fopen(tempFilename_.c_str(), "r");
    if (!in) {
      std::cerr << "Failed to open " << tempFilename_ << std::endl;
      return false;
    }

    std::FILE *out = fopen(newFilename.c_str(), "w");
    if (!out) {
      std::cerr << "Failed to init " << newFilename << std::endl;
      return false;
    }

    // Read one octree level at a time, starting from the back of source
    // Write to out file in reverse order

    // Position at end of file
    for (unsigned int ts=0; ts<numTimesteps_; ++ts) {

      size_t octreePos=static_cast<size_t>((numOTNodes)*numBrickVals*(ts+1));
      for (unsigned int level=0; level<numLevels_; ++level) {

        unsigned int bricksPerLevel = pow(8, level);
        unsigned int valuesPerLevel = numBrickVals*bricksPerLevel;
        octreePos -= valuesPerLevel;
        std::vector<float> buffer(valuesPerLevel);

        fseeko(in, octreePos*(size_t)sizeof(float), SEEK_SET);
        size_t readSize = static_cast<size_t>(valuesPerLevel)*sizeof(float);
        fread(reinterpret_cast<void*>(&buffer[0]), readSize, 1, in);
        fwrite(reinterpret_cast<void*>(&buffer[0]), readSize, 1, out);
      }

    }

    fclose(in);
    fclose(out);

  } // scoping streams

  // Create one file for every level of the BST tree structure
  // by averaging the values in the one below.
  unsigned int numTimestepsInLevel = numTimesteps_;
  unsigned int numValsInOT = numBrickVals*numOTNodes;
  std::vector<float> inBuffer1(numValsInOT);
  std::vector<float> inBuffer2(numValsInOT);
  std::vector<float> outBuffer(numValsInOT);

  size_t OTBytes = static_cast<size_t>(numValsInOT * sizeof(float));
  std::string fromFilename = newFilename;
  std::string toFilename;

  do {

    std::stringstream ss;
    ss << BSTLevel - 1;
    std::cout << "Creating level " << BSTLevel << std::endl;
    toFilename = tempFilename_ + "." + ss.str() + ".tmp";

    // Init files

    std::FILE *in = fopen(fromFilename.c_str(), "r");
    if (!in) {
      std::cerr << "Failed to open " << fromFilename << std::endl;
      return false;
    }

    std::FILE *out = fopen(toFilename.c_str(), "w");
    if (!out) {
      std::cerr << "Failed to init " << toFilename << std::endl;
      return false;
    }

    fseeko(in, 0, SEEK_END);
    fseeko(in, 0, SEEK_SET);

    for (unsigned int ts=0; ts<numTimestepsInLevel; ts+=2) {
    
      // Read two octrees (two time steps)
      fread(reinterpret_cast<void*>(&inBuffer1[0]), OTBytes, 1, in);
      fread(reinterpret_cast<void*>(&inBuffer2[0]), OTBytes, 1, in);

      // Average time steps
      for (unsigned int i=0; i<outBuffer.size(); ++i) {
        outBuffer[i] = (inBuffer1[i] + inBuffer2[i]) / static_cast<float>(2);
      }

      // Write brick
      fwrite(reinterpret_cast<void*>(&outBuffer[0]), OTBytes, 1, out);

    }

    fromFilename = toFilename;

    fclose(in);
    fclose(out);
    if (BSTLevel > 0) {
      BSTLevel--;
    }
    numTimestepsInLevel /= 2;
                
  } while (BSTLevel != 0);

  std::FILE *out = fopen(outFilename_.c_str(), "w");
  if (!out) {
    std::cerr << "Failed to open for write: " << outFilename_ << std::endl;
    return false;
  }

  // Write header
  std::cout << "Writing header" << std::endl;
  size_t s = sizeof(unsigned int);
  fwrite(reinterpret_cast<void*>(&gridType_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&origNumTimesteps), s, 1, out);
  fwrite(reinterpret_cast<void*>(&numTimesteps_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&xBrickDim_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&yBrickDim_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&zBrickDim_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&xNumBricks_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&yNumBricks_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&zNumBricks_), s, 1, out);

  for (unsigned int level=0; level<numBSTLevels; ++level) {
  
    std::cout << "Writing level " << level+1 << " to output" << std::endl;

    std::stringstream ss;
    ss << level;
    std::string fromFilename = tempFilename_ + "." + ss.str() + ".tmp";

    std::FILE *in = fopen(fromFilename.c_str(), "r");
    if (!in) {
      std::cerr << "Failed to open " << fromFilename << std::endl;
      return false;
    }
    
    fseeko(in, 0, SEEK_END);
    size_t inFileSize = ftello(in);
    fseeko(in, 0, SEEK_SET);
    
    size_t chunkSize = 1024 * 1024 * 512; // Write a maximum of 512 MB at a time
    std::vector<float> buffer((size_t)chunkSize/sizeof(float));

    for (size_t offset = 0; offset < inFileSize; offset += chunkSize) {
      size_t thisSize = std::min(chunkSize, inFileSize - offset);

      fread(reinterpret_cast<void*>(&buffer[0]),
              static_cast<size_t>(thisSize), 1, in);

      fwrite(reinterpret_cast<void*>(&buffer[0]),
              static_cast<size_t>(thisSize), 1, out);
    }

    fclose(in);
  }

  fclose(out);

  // Do some validation on the out file
  std::FILE *in = fopen(outFilename_.c_str(), "r");
  if (!in) {
    std::cerr << "Failed to open " << outFilename_ << " for validation" 
      << std::endl;
      return false;
  } 

  unsigned int gridTypeTest, origNumTimestepsTest, 
    numTimestepsTest, xBrickDimTest,
    yBrickDimTest, zBrickDimTest, xNumBricksTest, yNumBricksTest,
    zNumBricksTest; 
  s = sizeof(unsigned int);
  fread(reinterpret_cast<void*>(&gridTypeTest), s, 1, in);
  fread(reinterpret_cast<void*>(&origNumTimestepsTest), s, 1, in);
  fread(reinterpret_cast<void*>(&numTimestepsTest), s, 1, in);
  fread(reinterpret_cast<void*>(&xBrickDimTest), s, 1, in);
  fread(reinterpret_cast<void*>(&yBrickDimTest), s, 1, in);
  fread(reinterpret_cast<void*>(&zBrickDimTest), s, 1, in);
  fread(reinterpret_cast<void*>(&xNumBricksTest), s, 1, in);
  fread(reinterpret_cast<void*>(&yNumBricksTest), s, 1, in);
  fread(reinterpret_cast<void*>(&zNumBricksTest), s, 1, in);

  if (gridTypeTest != gridType_ ||
      origNumTimestepsTest != origNumTimesteps ||
      numTimestepsTest != numTimesteps_ ||
      xBrickDimTest != xBrickDim_ ||
      yBrickDimTest != yBrickDim_ ||
      zBrickDimTest != zBrickDim_ ||
      xNumBricksTest != xNumBricks_ ||
      yNumBricksTest != yNumBricks_ ||
      zNumBricksTest != zNumBricks_) {
    std::cerr << "\n\nWARNING: Header check failed!\n" << std::endl;
    std::cout << "Values from file: " << std::endl;
    std::cout << "Grid type: " << gridTypeTest << std::endl;
    std::cout << "Orig num of timesteps " << origNumTimestepsTest<<std::endl;
    std::cout << "Number of timesteps: " << numTimestepsTest << std::endl;
    std::cout << "Brick dimensions: " << xBrickDimTest << " x " 
                << yBrickDimTest << " x " << zBrickDimTest << std::endl;
    std::cout << "Number of bricks: " << xNumBricksTest << " x " 
              << yNumBricksTest << " x " << zNumBricksTest << std::endl;
  } else {
    std::cout << "Header OK!" << std::endl;
  }

  size_t dataPos = ftello(in);
  // Check file size
  fseek(in, 0, SEEK_END);
  size_t fileSize = ftello(in);
  size_t calcSize = static_cast<size_t>(numBricksTotal_) *
                 static_cast<size_t>(xPaddedBrickDim_*
                                  yPaddedBrickDim_*
                                  zPaddedBrickDim_) *
                 static_cast<size_t>(sizeof(float)) + dataPos;
  if (fileSize != calcSize) {
    std::cerr << "File sizes don't match" << std::endl;
    std::cerr << "Calculated file size: " << calcSize << std::endl;
    std::cerr << "Real file size: " << fileSize << std::endl;
    return false;
  } else {
    std::cout << "File sizes OK!" << std::endl;
  }
  fseek(in, dataPos, SEEK_SET);

  // TODO check for inf/NaN

  fclose(in);

  return true;
}

} // namespace osp


