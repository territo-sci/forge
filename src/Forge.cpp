/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#include <Forge.h>
#include <Brick.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <array>
#include <boost/filesystem.hpp>

using namespace osp;

Forge * Forge::New() {
  return new Forge();
}

Forge::Forge() 
  : inFilename_("NotSet"), outFilename_("NotSet"),
    structure_(0), brickDim_(0), paddingWidth_(0) {
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

void Forge::SetStructure(unsigned int _structure) {
  structure_ = _structure;
}

void Forge::SetBrickDimensions(unsigned int _brickDim) {
  brickDim_ = _brickDim;
}

void Forge::SetPaddingWidth(unsigned int _paddingWidth) {
  paddingWidth_ = _paddingWidth;
}

bool Forge::Construct() {

  // Read metadata 
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

  if (!ConstructTSPTreeSpatial()) {
    std::cout << "Error: Failed to construct TSP tree" << std::endl;
    return false;
  }


  /*
  // Use the octrees to construct a TSP tree, write directly to disk
  if (!ConstructTSPTree()) {
    std::cout << "Error: Failed to construct TSP tree" << std::endl;
    return false;
  }
  */
  
  // Delete temporary file
  //DeleteTempFile();

  return true;
}


bool Forge::ReadMetadata() {

  std::cout << "Reading metadata" << std::endl;

  if (brickDim_ == 0) {
    std::cout << "Warning: Brick dimensions are zero!"<< std::endl;
  }

  dataSize_ = sizeof(real);

  // Read from file
  if (instream_.is_open()) {
    std::cout << "Error: Instream is already open!" << std::endl;
    return false;
  }
  instream_.open(inFilename_.c_str(), 
                std::ios_base::in | std::ios_base::binary);
  if (!instream_.is_open()) {
    std::cout << "Error: Could not open file." << std::endl;
    return false;
  }
  
  size_t s = sizeof(unsigned int);
  instream_.read(reinterpret_cast<char*>(&dataDimensionality_), s);
  instream_.read(reinterpret_cast<char*>(&numTimesteps_), s);
  instream_.read(reinterpret_cast<char*>(&dim_), s);
  instream_.read(reinterpret_cast<char*>(&dim_), s);
  instream_.read(reinterpret_cast<char*>(&dim_), s);

  // TODO support non-full BST trees? Right now, the number of timesteps
  // needs to be a power of two. Abort if it's not.
  if ( (numTimesteps_ & (numTimesteps_-1)) != 0) {
    std::cout << "ERROR: Number of timesteps not  power of two" << std::endl;
    return false;
  }

  if (dim_ % brickDim_ !=0) {
    std::cout << "Error: Voxel and brick dimension mismatch!" << std::endl;
    instream_.close();
    return false;
  }

  numBricks_ = dim_ / brickDim_;
  paddedDim_ = dim_ + paddingWidth_*2;
  paddedBrickDim_ = brickDim_ + paddingWidth_*2;

  std::cout << std::endl << "FORGE METADATA" << std::endl;
  std::cout << "Read from " << inFilename_ << " complete!" << std::endl;
  std::cout << "Structure: " << structure_ << std::endl;
  std::cout << "Data dimensionality: " << dataDimensionality_ << std::endl;
  std::cout << "Number of timesteps: " << numTimesteps_ << std::endl;
  std::cout << "Dimensions: " << dim_ << " x " << dim_ <<  
               " x " << dim_ << std::endl;
  std::cout << "Brick dimensions: " << brickDim_ << " x " << brickDim_ <<
               " x " << brickDim_ << std::endl;
  std::cout << "Padded dimensions: " << paddedDim_ << " x " << paddedDim_ <<
               " x " << paddedDim_ << std::endl;
  std::cout << "Padded brick dimensions: " << paddedBrickDim_ << " x " <<
                paddedBrickDim_ << " x " << paddedBrickDim_ << std::endl;
  std::cout << "Number of bricks: " << numBricks_ << " x " << numBricks_ << 
               " x " << numBricks_ << std::endl;
  std::cout << "Data size (bytes): " << dataSize_ << std::endl;
  std::cout << "Out file name: " << outFilename_ << std::endl;

  // Calculate some common things

  // Number of bricks in the base (leaf) level
  nrBricksBaseLevel_ = numBricks_*numBricks_*numBricks_;
  // Number of octree levels
  nrLevels_ = log(numBricks_)/log(2) + 1;
  // Number of bricks per octree (used for offsets)
  nrBricksPerOctree_ = (pow(8, nrLevels_) - 1) / 7;

  std::cout << "Number of bricks in base octree level: " 
    << nrBricksBaseLevel_ << std::endl;
  std::cout << "Number of levels in octree: " << nrLevels_ << std::endl;
  std::cout << "Number of bricks in octree: " 
    << nrBricksPerOctree_ << std::endl;

  // Save position of first data entry after header
  headerOffset_ = instream_.tellg();
  
  instream_.close();

  return true;
}

bool Forge::CreateOctree() {
  
  // Init out file
  std::fstream out;
  out.open(tempFilename_.c_str(),
           std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
  if (!out.is_open()) {
    std::cout << "Could not open " << tempFilename_ << " for write"<<std::endl;
    return false;
  }
  unsigned int outpos = 0;

  // Init in file
  std::fstream in;
  in.open(inFilename_.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!in.is_open()) {
    std::cout << "Could not open " << inFilename_ <<" for reading"<<std::endl;
  }

  // Loop over all timesteps 
  for (unsigned int i=0; i<numTimesteps_; ++i) {
    
    std::cout << "Constructing octree for timestep " << i << std::endl;
    // Read whole timestep into memory

    std::vector<real> timestepData(dim_*dim_*dim_, static_cast<real>(0));

    // Point to the right position in the file stream
    unsigned int timestepSize = dim_*dim_*dim_*dataSize_;
    std::ios::pos_type timestepOffset = 
      static_cast<std::ios::pos_type>(i*timestepSize) + headerOffset_;
    in.seekg(timestepOffset);
    // Read data
    in.read(reinterpret_cast<char*>(&timestepData[0]), timestepSize);

    // We now have a non-padded time step, and need to pad the borders.
    // If the padding width is 0, this should result in the same volume as 
    // before.
    std::vector<real> paddedData(paddedDim_*paddedDim_*paddedDim_,
                                 static_cast<real>(0));
    // Loop over padded volume to fill
    // xp -> "x padded"
    // xo -> "x original"
    unsigned int xo, yo, zo;
    for (unsigned int zp=0; zp<paddedDim_; ++zp) {
      for (unsigned int yp=0; yp<paddedDim_; ++yp) {
        for (unsigned int xp=0; xp<paddedDim_; ++xp) {
          
          if (xp == 0) {
            xo = xp;
          } else if (xp == paddedDim_-1) {
            xo = xp-2;
          } else {
            xo = xp-1;
          }

          if (yp == 0) {
            yo = yp;
          } else if (yp == paddedDim_-1) {
            yo = yp-2;
          } else {
            yo = yp-1;
          }
          
          if (zp == 0) {
            zo = zp;
          } else if (zp == paddedDim_-1) {
            zo = zp-2;
          } else {
            zo = zp-1;
          }

          paddedData[xp + yp*paddedDim_ + zp*paddedDim_*paddedDim_] =
            timestepData[xo + yo*dim_ + zo*dim_*dim_];

        }
      }
    }

    // Create a container for the base level bricks.
    std::vector<Brick<real>* > baseLevelBricks(nrBricksBaseLevel_, NULL);

    // Loop over the volume's subvolumes and create one brick for each
    for (unsigned int zBrick=0; zBrick<numBricks_; ++zBrick) {
      for (unsigned int yBrick=0; yBrick<numBricks_; ++yBrick) {
        for (unsigned int xBrick=0; xBrick<numBricks_; ++xBrick) {
        
          Brick<real> *brick = Brick<real>::New(paddedBrickDim_, 
                                                paddedBrickDim_, 
                                                paddedBrickDim_,
                                                static_cast<real>(0));  
            
          // Loop over the subvolume's voxels
          unsigned int xMin = xBrick * brickDim_;
          unsigned int xMax = (xBrick + 1) * brickDim_ - 1 + paddingWidth_*2;
          unsigned int yMin = yBrick * brickDim_;
          unsigned int yMax = (yBrick + 1) * brickDim_ - 1 + paddingWidth_*2;
          unsigned int zMin = zBrick * brickDim_;
          unsigned int zMax = (zBrick + 1) * brickDim_ - 1 + paddingWidth_*2;
          //std::cout << "xmin / xmax " << xMin << " / " << xMax << std::endl;
          unsigned int zLoc= 0;
          for (unsigned int zSub=zMin; zSub<=zMax; ++zSub) {
            unsigned int yLoc = 0;
            for (unsigned int ySub=yMin; ySub<=yMax; ++ySub) {  
              unsigned int xLoc = 0;
              for (unsigned int xSub=xMin; xSub<=xMax; ++xSub) {
                // Look up global index in full volume
                unsigned int globalIndex = 
                  xSub + ySub*paddedDim_ + zSub*paddedDim_*paddedDim_;
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
            xBrick + yBrick*numBricks_ + zBrick*numBricks_*numBricks_;
          baseLevelBricks[brickIndex] = brick;
        }
      }
    }

    // Make a container for all the octree bricks
    std::vector<Brick<real>* > octreeBricks(nrBricksPerOctree_);
    // Use Z-order coordinates to rearrange the base level bricks

    // so that the eight children for each parent node lies
    // next to each other
    for (uint16_t z=0; z<(uint16_t)numBricks_; ++z) { 
      for (uint16_t y=0; y<(uint16_t)numBricks_; ++y) {
        for (uint16_t x=0; x<(uint16_t)numBricks_; ++x) {
          unsigned int zOrderIdx = static_cast<unsigned int>(ZOrder(x, y, z));
          unsigned int idx = x + y*numBricks_ + z*numBricks_*numBricks_;
          octreeBricks[zOrderIdx] = baseLevelBricks[idx];

        }
      }
    }

    // Construct higher levels of octree

    // Position for next brick, starting at position beyond base level
    unsigned int brickPos = nrBricksBaseLevel_;
    // Position for first child to average
    unsigned int childPos = 0;

    while (brickPos < nrBricksPerOctree_) {

      // Filter the eight children and then combine them to build
      // the higher level node
      std::vector<Brick<real>* > filteredChildren(8, NULL);
      unsigned int i=0;
      for (unsigned int child=childPos; child<childPos+8; ++child) {
        Brick<real> *filteredChild = 
          Brick<real>::Filter(octreeBricks[child]);
        filteredChildren[i++] = filteredChild;
      }
      Brick<real> *newBrick = Brick<real>::Combine(filteredChildren);
      
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
    int n = 0;
    for (auto it=octreeBricks.begin(); it!=octreeBricks.end(); ++it) {
      out.write(reinterpret_cast<char*>(&(*it)->data_[0]), (*it)->Size());
      // Free memory when we're done
      delete *it;
    }
    
  }  
 
  in.close();
  out.close();

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


bool Forge::DeleteTempFile() {
  if (boost::filesystem::exists(tempFilename_)) {
    boost::filesystem::remove_all(tempFilename_);
  } else {
    std::cout <<"Warning: " << tempFilename_ << " does not exist" << std::endl;
  }
  return true;
}

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
  unsigned int numOctreeNodes = nrBricksPerOctree_;
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
  while (level < nrLevels_) {

    unsigned int bricksPerLevel = pow(8, level);
    unsigned int valuesPerLevel = brickSize * bricksPerLevel;
    octreePos -= valuesPerLevel;
    //std::cout << "Level " << level << 
    //  ", starting octree pos: " << octreePos << std::endl;
    
    // Loop over all octree nodes in this level
    for (unsigned int i=0; i<bricksPerLevel; ++i) {

      // Allocate bricks for one BST
      std::vector<Brick<real>* > BSTBricks(numBSTNodes);

      // Starting position in the BST (base level first)
      // Note that we want the root to be in the front
      unsigned int BSTBrickPos = numBSTNodes - numTimesteps_;

      // Collect all corresponding bricks from the octree and build the leaves
      for (unsigned int ts=0; ts<numTimesteps_; ++ts) {
        Brick<real> *brick = Brick<real>::New(paddedBrickDim_,
                                              paddedBrickDim_,
                                              paddedBrickDim_,
                                              static_cast<real>(0));
        real *dataPtr = &(brick->data_[0]);
        size_t s = brickSize*sizeof(real);
        in.seekg(static_cast<std::ios::pos_type>(
          (octreePos+ts*numOctreeNodes*brickSize)*sizeof(real)));
        in.read(reinterpret_cast<char*>(dataPtr), s);

        BSTBricks[BSTBrickPos++] = brick;
      }
      // Rewind position
      BSTBrickPos -= numTimesteps_;

      // Average bricks to build higher levels in BST
      // This is really a reversed level, but it works well for loop purposes
      unsigned int level = 1;
      do {
        // Save position for bricks to average from
        unsigned int fromPos = BSTBrickPos;
        // Calculate starting position
        unsigned int bricksAtLevel = numTimesteps_/(pow(2, level));
        BSTBrickPos -= bricksAtLevel;
        // Average bricks 
        for (int i=0; i<bricksAtLevel; ++i) {
          BSTBricks[BSTBrickPos] = Brick<real>::Average(BSTBricks[fromPos],
                                                        BSTBricks[fromPos+1]);
          BSTBrickPos++;
          fromPos += 2;
        }
        level++;
        BSTBrickPos -= bricksAtLevel;
      } while (BSTBrickPos > 0);
      
      // Write BST to file
      unsigned int s = brickSize*sizeof(real);
      for (auto it=BSTBricks.begin(); it!=BSTBricks.end(); ++it) {
        out.write(reinterpret_cast<char*>(&(*it)->data_[0]), s);   
      }

      octreePos += brickSize;
 
    }
    // Rewind
    octreePos -= valuesPerLevel;

    level++;

  } // while level < nrLevels

  return true;
} 


bool Forge::ConstructTSPTreeSpatial() {

  // Make sure temporary file exists
  if (!boost::filesystem::exists(tempFilename_)) {
    std::cout << "Error: temp file "<<tempFilename_<<" missing" << std::endl;
    return false;
  }

  // Numbers to keep track of
  unsigned int numOTNodes = nrBricksPerOctree_;
  unsigned int numBSTNodes = 2*numTimesteps_ - 1;
  unsigned int numBSTLevels = log(numTimesteps_)/log(2) + 1;
  unsigned int numBrickVals = 
    paddedBrickDim_*paddedBrickDim_*paddedBrickDim_;

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
    std::cout << "Creating " << newFilename << std::endl;
    boost::filesystem::rename(tempFilename_, newFilename);
  }

  // Create one file for every level of the BST tree structure
  // by averaging the values in the one below.

  // Read two bricks at a time and average
  std::vector<real> inBuffer1(numBrickVals);
  std::vector<real> inBuffer2(numBrickVals);
  std::vector<real> outBuffer(numBrickVals);
  unsigned int brickSize = numBrickVals * sizeof(real);

  unsigned int numTimestepsInLevel = numTimesteps_;
  std::string fromFilename = newFilename;
  std::string toFilename;
  do {

    std::stringstream ss;
    ss << BSTLevel - 1;
    toFilename = tempFilename_ + "." + ss.str() + ".tmp";
    
    // Init in stream
    std::fstream in;
    in.open(fromFilename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!in.is_open()) {
      std::cout << "Error: Could not open " << fromFilename << std::endl;
      return false;
    }
    
    // Init out stream
    std::fstream out;
    out.open(toFilename.c_str(),
            std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    if (!out.is_open()) {
      std::cout << "Could not open " << toFilename << " for write"<<std::endl;
      return false;
    }
    std::cout << "Writing to " << toFilename << std::endl;

    for (unsigned int ts=0; ts<numTimestepsInLevel; ts+=2) {
    
      unsigned int fromPos1 = ts*numBrickVals;
      unsigned int fromPos2 = (ts+1)*numBrickVals;
      unsigned int toPos = ts/2;

      for (unsigned int brick=0; brick<numOTNodes; ++brick) {

        // Read two time steps
        in.read(reinterpret_cast<char*>(&inBuffer1[0]), brickSize);
        in.read(reinterpret_cast<char*>(&inBuffer2[0]), brickSize);
        
        // Average time steps
        for (unsigned int i=0; i<outBuffer.size(); ++i) {
          outBuffer[i] = (inBuffer1[i] + inBuffer2[i]) / static_cast<real>(2);
        }

        // Write brick
        out.write(reinterpret_cast<char*>(&outBuffer[0]), brickSize);

      }

    }

    fromFilename = toFilename;
    in.close();
    out.close();

    BSTLevel--;
    numTimestepsInLevel /= 2;
                
  } while (BSTLevel != 0);

  // The octrees are now in the temporary files, and all we need to do
  // is reading them back the right order and write to the output

  // Open out filestream
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
  
  for (unsigned int level=0; level<numBSTLevels; ++level) {
  
    std::cout << "Writing level " << level << " to output" << std::endl;

    std::stringstream ss;
    ss << level;
    std::string fromFilename = tempFilename_ + "." + ss.str() + ".tmp";

    std::fstream in;
    in.open(fromFilename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!in.is_open()) {
      std::cout << "Error: Could not open " << fromFilename << std::endl;
      return false;
    }
  
    // Check file size correctness
    in.seekg(0, std::ifstream::end);
    unsigned int realSize = in.tellg();
    unsigned int numTimesteps = pow(2, level);
    unsigned int expectedSize = 
      numBrickVals*numOTNodes*numTimesteps*sizeof(real);
    if (realSize != expectedSize) {
      std::cout << "Error: File size not as expected" << std::endl;
      return false;
    }

    // Read whole file, write to out file
    std::vector<real> buffer(realSize/sizeof(real));
    in.read(reinterpret_cast<char*>(&buffer[0]), realSize);
    out.write(reinterpret_cast<char*>(&buffer[0]), realSize);

    in.close();
  }

  out.close();

  return true;
}
