/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 * Note: Uses C-style I/O to handle large files. This functionality
 * might be implementation/platform dependent!
 *
 */

#include <Forge.h>
#include <Brick.h>
#include <iostream>
#include <sstream>
//#include <fstream>
#include <math.h>
#include <array>
#include <boost/filesystem.hpp>
#include <stdio.h>

using namespace osp;

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
    zPaddedBrickDim_(0),
    spatialScaling_(1.f),
    temporalScaling_(1.f) {
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

void Forge::SetSpatialScaling(float _spatialScaling) {
  spatialScaling_ = _spatialScaling;
}

void Forge::SetTemporalScaling(float _temporalScaling) {
  temporalScaling_ = _temporalScaling;
}

bool Forge::Construct() {

  // Read metadata from VDF file
  std::cout << "Reading metadata" << std::endl;
  if (!ReadMetadata()) {
    std::cout << "Error: Could not read metadata" << std::endl;
    return false;
  }

  // Create one octree per timestep, save in temporary file
  if (!CreateOctree()) {
    std::cout << "Error: Failed to create temp octree file" << std::endl;
    return false;
  }

  // Construct the TSP tree from the temporary files
  if (!ConstructTSPTree()) {
    std::cout << "Error: Failed to construct TSP tree" << std::endl;
    return false;
  }

  // Delete temporary files
  if (!DeleteTempFiles()) {
    std::cout << "Failed to delete temp files, but that's okay" << std::endl;
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
  if ( (numTimesteps_ & (numTimesteps_-1)) != 0) {
    std::cerr << "ERROR: Number of timesteps not power of two" << std::endl;
    fclose(in);
    return false;
  }

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
  headerOffset_ = ftell(in);
  
  fclose(in);

  return true;
}


bool Forge::CreateOctree() {
  
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
  for (unsigned int i=0; i<numTimesteps_; ++i) {
    
    std::cout << "Constructing octree for timestep " << i << "/" 
      << numTimesteps_ << "\r" << std::flush;

    std::vector<float> timestepData(xDim_*yDim_*zDim_, static_cast<float>(0));

    // Point to the right position in the file stream and read it
    size_t timestepSize = xDim_*yDim_*zDim_*dataSize_;
    size_t timestepOffset = static_cast<size_t>(i)*timestepSize+headerOffset_;
    fseek(in, timestepOffset, SEEK_SET);
    fread(reinterpret_cast<void*>(&timestepData[0]), timestepSize, 1, in);

    // We now have a non-padded time step, and need to pad the borders.
    // If the padding width is 0, this should result in the same volume as 
    // before.
    std::vector<float> paddedData(xPaddedDim_*yPaddedDim_*zPaddedDim_,
                                 static_cast<float>(0));
    // Loop over padded volume that we want to fill
    // xp -> "x padded"
    // xo -> "x original"
    unsigned int xo, yo, zo;
    for (unsigned int zp=0; zp<zPaddedDim_; ++zp) {
      for (unsigned int yp=0; yp<yPaddedDim_; ++yp) {
        for (unsigned int xp=0; xp<xPaddedDim_; ++xp) {
          
          if (xp == 0) {
            xo = xp;
          } else if (xp == xPaddedDim_-1) {
            xo = xp-2;
          } else {
            xo = xp-1;
          }

          if (yp == 0) {
            yo = yp;
          } else if (yp == yPaddedDim_-1) {
            yo = yp-2;
          } else {
            yo = yp-1;
          }
          
          if (zp == 0) {
            zo = zp;
          } else if (zp == zPaddedDim_-1) {
            zo = zp-2;
          } else {
            zo = zp-1;
          }

          paddedData[xp + yp*xPaddedDim_ + zp*xPaddedDim_*yPaddedDim_] =
            timestepData[xo + yo*xDim_ + zo*yDim_*zDim_];
        }
      }
    }

    std::cout<<"                                             \r"<<std::flush;

    // Create a container for the base level bricks.
    std::vector<Brick<float>* > baseLevelBricks(numBricksBaseLevel_, NULL);

    // Loop over the volume's subvolumes and create one brick for each
    for (unsigned int zBrick=0; zBrick<zNumBricks_; ++zBrick) {
      for (unsigned int yBrick=0; yBrick<yNumBricks_; ++yBrick) {
        for (unsigned int xBrick=0; xBrick<xNumBricks_; ++xBrick) {
        
          Brick<float> *brick = Brick<float>::New(xPaddedBrickDim_, 
                                                  yPaddedBrickDim_, 
                                                  zPaddedBrickDim_,
                                                  static_cast<float>(0));  
            
          // Loop over the subvolume's voxels
          unsigned int xMin = xBrick * xBrickDim_;
          unsigned int xMax = (xBrick + 1) * xBrickDim_ - 1 + paddingWidth_*2;
          unsigned int yMin = yBrick * yBrickDim_;
          unsigned int yMax = (yBrick + 1) * yBrickDim_ - 1 + paddingWidth_*2;
          unsigned int zMin = zBrick * zBrickDim_;
          unsigned int zMax = (zBrick + 1) * zBrickDim_ - 1 + paddingWidth_*2;
          //std::cout << "xmin / xmax " << xMin << " / " << xMax << std::endl;
          unsigned int zLoc= 0;
          for (unsigned int zSub=zMin; zSub<=zMax; ++zSub) {
            unsigned int yLoc = 0;
            for (unsigned int ySub=yMin; ySub<=yMax; ++ySub) {  
              unsigned int xLoc = 0;
              for (unsigned int xSub=xMin; xSub<=xMax; ++xSub) {
                // Look up global index in full volume
                unsigned int globalIndex = 
                  xSub + ySub*xPaddedDim_ + zSub*xPaddedDim_*yPaddedDim_;
                // Set data at local subvolume index
                brick->SetData(xLoc, yLoc, zLoc, paddedData[globalIndex]); 
                xLoc++;
              }
              yLoc++;
            }
            zLoc++;
          }

          // Save to base level
          unsigned int brickIndex = 
            xBrick + yBrick*xNumBricks_ + zBrick*xNumBricks_*yNumBricks_;
          baseLevelBricks[brickIndex] = brick;
        }
      }
    }

    // Make a container for all the octree bricks
    std::vector<Brick<float>* > octreeBricks(numBricksPerOctree_);
    // Use Z-order coordinates to rearrange the base level bricks

    // so that the eight children for each parent node lies
    // next to each other
    for (uint16_t z=0; z<static_cast<uint16_t>(xNumBricks_); ++z) { 
      for (uint16_t y=0; y<static_cast<uint16_t>(yNumBricks_); ++y) {
        for (uint16_t x=0; x<static_cast<uint16_t>(zNumBricks_); ++x) {
          unsigned int zOrderIdx = static_cast<unsigned int>(ZOrder(x, y, z));
          unsigned int idx = x + y*xNumBricks_ + z*xNumBricks_*yNumBricks_;
          octreeBricks[zOrderIdx] = baseLevelBricks[idx];
        }
      }
    }

    // Construct higher levels of octree

    // Position for next brick, starting at position beyond base level
    unsigned int brickPos = numBricksBaseLevel_;
    // Position for first child to average
    unsigned int childPos = 0;


    while (brickPos < numBricksPerOctree_) {
      // Filter the eight children and then combine them to build
      // the higher level node
      std::vector<Brick<float>* > filteredChildren(8, NULL);
      unsigned int i=0;
      for (unsigned int child=childPos; child<childPos+8; ++child) {
        Brick<float> *filteredChild = 
          Brick<float>::Filter(octreeBricks[child]);
        filteredChildren[i++] = filteredChild;
      }
      Brick<float> *newBrick = Brick<float>::Combine(filteredChildren);
      
      for (auto it=filteredChildren.begin(); it!=filteredChildren.end(); ++it){
        delete *it;
        *it = NULL;
      }

      // Set next child pos
      childPos += 8;

      // Save new brick
      octreeBricks[brickPos++] = newBrick;
    }

    // Write octree to file
    for (auto it=octreeBricks.begin(); it!=octreeBricks.end(); ++it) {
      fwrite(reinterpret_cast<void*>(&(*it)->data_[0]),
             static_cast<size_t>((*it)->Size()), 1, out);
      // Free memory when we're done
      delete *it;
    }
    
  }  
 
  fclose(in);
  fclose(out);

  return true;
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

/*
bool Forge::ConstructTSPTree() {

  // Make sure the temporary file exists
  if (!boost::filesystem::exists(tempFilename_)) {
    std::cout << "Error: Temp file "<<tempFilename_ <<" missing"<< std::endl;
    return false;
  }
  std::fstream in;
  in.open(tempFilename_.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!in.is_open()) {
    std::cout << "Error: Could not open " << tempFilename_ << std::endl;
    return false;
  }

  // Number of nodes in the octree skeleton
  unsigned int numOctreeNodes = numBricksPerOctree_;
  std::cout << "Number of nodes per octree " << numOctreeNodes << std::endl;

  // Number of nodes in binary time tree
  unsigned int numBSTNodes = 2*numTimesteps_ - 1;
  std::cout << "Num nodes in BST: " << numBSTNodes << std::endl;

  std::fstream out;
  out.open(tspFilename_.c_str(), 
           std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
  if (!out.is_open()) {
    std::cout << "Error: could not open " << tspFilename_ << std::endl;
    return false;
  }
  
  // Write header
  out.seekp(std::ios_base::beg);
  std::cout << "Writing header" << std::endl;
  // Write header
  size_t s = sizeof(unsigned int);
  out.write(reinterpret_cast<char*>(&structure_), s);
  out.write(reinterpret_cast<char*>(&dataDimensionality_), s);
  out.write(reinterpret_cast<char*>(&brickDim_), s);
  out.write(reinterpret_cast<char*>(&brickDim_), s);
  out.write(reinterpret_cast<char*>(&brickDim_), s);
  out.write(reinterpret_cast<char*>(&numBricks_), s);
  out.write(reinterpret_cast<char*>(&numBricks_), s);
  out.write(reinterpret_cast<char*>(&numBricks_), s);
  out.write(reinterpret_cast<char*>(&numTimesteps_), s);
  out.write(reinterpret_cast<char*>(&paddingWidth_), s);
  out.write(reinterpret_cast<char*>(&dataSize_), s);

  std::cout << "Position after writing header: " << out.tellp() << std::endl;

  // Loop over all positions in octree skeleton
  
  // Octree level (to both read from and write to)
  unsigned int level = 0;

  // Number of values in a brick
  unsigned int brickSize = paddedBrickDim_*paddedBrickDim_*paddedBrickDim_; 
  std::cout << "Num values in brick " << brickSize << std::endl;

  // Position in octree to read from (starting at first timestep)
  // Since the octree is constructed from the base level up,
  // the root is at the last brick index
  unsigned int octreePos = (numOctreeNodes)*brickSize;

  // Loop over all octree levels (in backwards order)
  while (level < numLevels_) {

    unsigned int bricksPerLevel = pow(8, level);
    unsigned int valuesPerLevel = brickSize * bricksPerLevel;
    octreePos -= valuesPerLevel;
    //std::cout << "Level " << level << 
    //  ", starting octree pos: " << octreePos << std::endl;
    
  unsigned int numBSTLevels = log(numTimesteps_)/log(2) + 1;
    // Loop over all octree nodes in this level
    for (unsigned int i=0; i<bricksPerLevel; ++i) {

      // Allocate bricks for one BST
      std::vector<Brick<float>* > BSTBricks(numBSTNodes);

      // Starting position in the BST (base level first)
      // Note that we want the root to be in the front
      unsigned int BSTBrickPos = numBSTNodes - numTimesteps_;

      // Collect all corresponding bricks from the octree and build the leaves
      for (unsigned int ts=0; ts<numTimesteps_; ++ts) {
        Brick<float> *brick = Brick<float>::New(paddedBrickDim_,
                                              paddedBrickDim_,
                                              paddedBrickDim_,
                                              static_cast<float>(0));
        float *dataPtr = &(brick->data_[0]);
        size_t s = brickSize*sizeof(float);
        in.seekg(static_cast<std::ios::pos_type>(
          (octreePos+ts*numOctreeNodes*brickSize)*sizeof(float)));
        in.read(reinterpret_cast<char*>(dataPtr), s);

        BSTBricks[BSTBrickPos++] = brick;
      }
      // Rewind position
      BSTBrickPos -= numTimesteps_;

      // Average bricks to build higher levels in BST
      // This is floatly a reversed level, but it works well for loop purposes
      unsigned int level = 1;
      do {
        // Save position for bricks to average from
        unsigned int fromPos = BSTBrickPos;
        // Calculate starting position
        unsigned int bricksAtLevel = numTimesteps_/(pow(2, level));
        BSTBrickPos -= bricksAtLevel;
        // Average bricks 
        for (int i=0; i<bricksAtLevel; ++i) {
          BSTBricks[BSTBrickPos] = Brick<float>::Average(BSTBricks[fromPos],
                                                        BSTBricks[fromPos+1]);
          BSTBrickPos++;
          fromPos += 2;
        }
        level++;
        BSTBrickPos -= bricksAtLevel;
      } while (BSTBrickPos > 0);
      
      // Write BST to file
      unsigned int s = brickSize*sizeof(float);
      for (auto it=BSTBricks.begin(); it!=BSTBricks.end(); ++it) {
        out.write(reinterpret_cast<char*>(&(*it)->data_[0]), s);   
      }

      octreePos += brickSize;
 
    }
    // Rewind
    octreePos -= valuesPerLevel;

    level++;

  } // while level < numLevels

  return true;
} 
*/

bool Forge::ConstructTSPTree() {

  // Make sure temporary file exists
  if (!boost::filesystem::exists(tempFilename_)) {
    std::cout << "Error: temp file "<<tempFilename_<<" missing" << std::endl;
    return false;
  }

  // Numbers to keep track of
  unsigned int numOTNodes = numBricksPerOctree_;
  unsigned int numBSTNodes = 2*numTimesteps_ - 1;
  unsigned int numBSTLevels = log(numTimesteps_)/log(2) + 1;
  unsigned int numBrickVals = 
    xPaddedBrickDim_*yPaddedBrickDim_*zPaddedBrickDim_;

  std::cout << "Num nodes per OT: " << numOTNodes << std::endl;
  std::cout << "Num nodes per BST: " << numBSTNodes << std::endl;
  std::cout << "Num BST levels: " << numBSTLevels << std::endl;
  std::cout << "Num values per brick: " << numBrickVals << std::endl;

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

      size_t octreePos = static_cast<size_t>((numOTNodes)*numBrickVals*(ts+1));
      for (unsigned int level=0; level<numLevels_; ++level) {

        unsigned int bricksPerLevel = pow(8, level);
        unsigned int valuesPerLevel = numBrickVals*bricksPerLevel;
        octreePos -= valuesPerLevel;
        std::vector<float> buffer(valuesPerLevel);

        fseek(in, octreePos*sizeof(float), SEEK_SET);
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

    //std::cout << "Writing to " << toFilename << std::endl;
    //std::cout << "From " << fromFilename << std::endl;

    fseek(in, 0, SEEK_END);
    size_t fileSize = ftell(in);
    //std::cout << "In file size: " << fileSize << std::endl;

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

    BSTLevel--;
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
  fwrite(reinterpret_cast<void*>(&numTimesteps_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&xBrickDim_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&yBrickDim_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&zBrickDim_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&xNumBricks_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&yNumBricks_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&zNumBricks_), s, 1, out);
  fwrite(reinterpret_cast<void*>(&dataSize_), s, 1, out);

  for (unsigned int level=0; level<numBSTLevels; ++level) {
  
    std::cout << "Writing level " << level+1 << " to output" << std::endl;

    std::stringstream ss;
    ss << level;
    std::string fromFilename = tempFilename_ + "." + ss.str() + ".tmp";
    //std::cout << "Reading from: " << fromFilename << std::endl;

    std::FILE *in = fopen(fromFilename.c_str(), "r");
    if (!in) {
      std::cerr << "Failed to open " << fromFilename << std::endl;
      return false;
    }
    
    fseek(in, 0, SEEK_END);
    size_t inFileSize = ftell(in);
    fseek(in, 0, SEEK_SET);
    
    std::vector<float> buffer(inFileSize/sizeof(float));
    // Read whole file, write to out file
    //in.read(reinterpret_cast<char*>(&buffer[0]), floatSize);
    fread(reinterpret_cast<void*>(&buffer[0]), inFileSize, 1, in);

    //out.write(reinterpret_cast<char*>(&buffer[0]), floatSize); 
    fwrite(reinterpret_cast<void*>(&buffer[0]), inFileSize, 1, out);
    //std::cout << "Pos after writing: " << ftell(out) << std::endl;

    fclose(in);
  }

  fclose(out);

  /*
  { // Scoping stream

    // Do some validation on the out file
    std::fstream test;
    test.open(outFilename_.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!test.is_open()) {
      std::cout << "Error: could not open test:" << outFilename_ << std::endl;
      return false;
    }
    
    // Read the header
    std::cout << "Validating header" << std::endl;
    unsigned int gridTypeTest, numTimestepsTest, xBrickDimTest,
      yBrickDimTest, zBrickDimTest, xNumBricksTest, yNumBricksTest,
      zNumBricksTest, dataSizeTest;
    test.read(reinterpret_cast<char*>(&gridTypeTest), s);
    test.read(reinterpret_cast<char*>(&numTimestepsTest), s);
    test.read(reinterpret_cast<char*>(&xBrickDimTest), s);
    test.read(reinterpret_cast<char*>(&yBrickDimTest), s);
    test.read(reinterpret_cast<char*>(&zBrickDimTest), s);
    test.read(reinterpret_cast<char*>(&xNumBricksTest), s);
    test.read(reinterpret_cast<char*>(&yNumBricksTest), s);
    test.read(reinterpret_cast<char*>(&zNumBricksTest), s);
    test.read(reinterpret_cast<char*>(&dataSizeTest), s);

    if (gridTypeTest != gridType_ ||
        numTimestepsTest != numTimesteps_ ||
        xBrickDimTest != xBrickDim_ ||
        yBrickDimTest != yBrickDim_ ||
        zBrickDimTest != zBrickDim_ ||
        xNumBricksTest != xNumBricks_ ||
        yNumBricksTest != yNumBricks_ ||
        zNumBricksTest != zNumBricks_ ||
        dataSizeTest != dataSizeTest) {
      std::cerr << "\n\nWARNING: Header check failed!\n" << std::endl;
      std::cout << "Values from file: " << std::endl;
      std::cout << "Grid type: " << gridTypeTest << std::endl;
      std::cout << "Number of timesteps: " << numTimestepsTest << std::endl;
      std::cout << "Brick dimensions: " << xBrickDimTest << " x " 
                  << yBrickDimTest << " x " << zBrickDimTest << std::endl;
      std::cout << "Number of bricks: " << xNumBricksTest << " x " 
                << yNumBricksTest << " x " << zNumBricksTest << std::endl;
      std::cout << "Data size (bytes): " << dataSizeTest << std::endl;
    } else {
      std::cout << "Header OK!" << std::endl;
    }

    // Read one brick at a time, check for NaN's or inf's
    std::cout << "Checking " << numBricksTotal_ << " bricks..." << std::endl;
    std::vector<float> brickBuffer(numBrickVals);
    for (unsigned int i=0; i<numBricksTotal_; ++i) {
      test.read(reinterpret_cast<char*>(&brickBuffer[0]), 
                  sizeof(float)*numBrickVals);
      for (auto it=brickBuffer.begin(); it!=brickBuffer.end(); ++it) {
        if (isnan(*it)) {
          std::cerr << "NaN detected: " << *it << std::endl;
        } 
        if (isinf(*it)) {
          std::cerr << "inf detected: " << *it << std::endl;
        }
      }
    }
    std::cout << "inf/NaN check complete" << std::endl;


    test.close();

  } // Scoping stream
  */
  


  return true;
}




