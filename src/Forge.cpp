/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#include <Forge.h>
#include <BricksHeader.h>
#include <Brick.h>
#include <iostream>
#include <fstream>

using namespace osp;

Forge * Forge::New() {
  return new Forge();
}

Forge::Forge() 
  : inFilename_("NotSet"), outFilename_("NotSet"), header_(NULL),
    structure_(0), xBrickDim_(0), yBrickDim_(0), zBrickDim_(0) {
}

Forge::~Forge() {
  if (header_) {
    delete header_;
  }
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

void Forge::SetBrickDimensions(unsigned int _xBrickDim,
                               unsigned int _yBrickDim,
                               unsigned int _zBrickDim) {
  xBrickDim_ = _xBrickDim;
  yBrickDim_ = _yBrickDim;
  zBrickDim_ = _zBrickDim;
}

bool Forge::Read() {

  // Read header information
  if (header_) {
    std::cout << "Warning: Header already exists, deleting it!" << std::endl;
    delete header_;
  }

  header_ = BricksHeader::New();
  
  // Set the data not read from file
  header_->SetStructure(structure_);
  if (xBrickDim_ == 0 || yBrickDim_ == 0 || zBrickDim_ == 0) {
    std::cout << "Warning: One or more brick dimensions are zero!"<< std::endl;
  }
  header_->SetDimensions(xBrickDim_, yBrickDim_, zBrickDim_);
  header_->SetDataSize(static_cast<unsigned int>(sizeof(real)));

  // Read from file
  std::fstream instream;
  instream.open(inFilename_.c_str(), 
                std::ios_base::in | std::ios_base::binary);
  if (!instream.is_open()) {
    std::cout << "Error: Could not open file." << std::endl;
    delete header_;
    return false;
  }
  
  unsigned int dataDimensionality = 0, numTimesteps = 0,
               xDim = 0, yDim = 0, zDim = 0;
  size_t s = sizeof(unsigned int);

  instream.read(reinterpret_cast<char*>(&dataDimensionality), s);
  instream.read(reinterpret_cast<char*>(&numTimesteps), s);
  instream.read(reinterpret_cast<char*>(&xDim), s);
  instream.read(reinterpret_cast<char*>(&yDim), s);
  instream.read(reinterpret_cast<char*>(&zDim), s);

  if (xDim % xBrickDim_ !=0 ||yDim % yBrickDim_ !=0||zDim % zBrickDim_ !=0) {
    std::cout << "Error: Voxel and brick dimension mismatch!" << std::endl;
    instream.close();
    delete header_;
    return false;
  }

  unsigned int xNumBricks = xDim / xBrickDim_;
  unsigned int yNumBricks = yDim / yBrickDim_;
  unsigned int zNumBricks = zDim / zBrickDim_;
  header_->SetNumBricks(xNumBricks, yNumBricks, zNumBricks);
  bricks_.resize(xNumBricks*yNumBricks*zNumBricks);

  std::cout << "Read from " << inFilename_ << " complete!" << std::endl;
  std::cout << "Data dimensionality: " << dataDimensionality << std::endl;
  std::cout << "Number of timesteps: " << numTimesteps << std::endl;
  std::cout << "Dimensions: " << xDim << " x " << yDim <<  
               " x " << zDim << std::endl;
  std::cout << "Brick dimensions: " << xBrickDim_ << " x " << yBrickDim_ <<
               " x " << zBrickDim_ << std::endl;
  std::cout << "Number of bricks: " << xNumBricks << " x " << yNumBricks << 
               " x " << zNumBricks << std::endl;
  std::cout << "Structure: " << structure_ << std::endl;
  std::cout << "Data size (bytes): " << header_->DataSize() << std::endl;
  std::cout << "Out file name: " << outFilename_ << std::endl;

  // Read bricks
  for (unsigned int timestep=0; timestep<numTimesteps; ++timestep) {
  
    // Read whole timestep into memory
    unsigned int size = xDim*yDim*zDim*sizeof(float);
    std::vector<real> timestepData(xDim*yDim*zDim, static_cast<real>(0));
    instream.read(reinterpret_cast<char*>(&timestepData[0]), size);

    // Loop over the volume's subvolumes and create one brick for each
    for (unsigned int zBrick=0; zBrick<zNumBricks; ++zBrick) {
      for (unsigned int yBrick=0; yBrick<yNumBricks; ++yBrick) {
        for (unsigned int xBrick=0; xBrick<xNumBricks; ++xBrick) {
        
          Brick<real> *brick = Brick<real>::New(xBrickDim_, 
                                                yBrickDim_, 
                                                zBrickDim_,
                                                static_cast<real>(0));  
            
          // Loop over the subvolume's voxels
          unsigned int xMin = xBrick * xBrickDim_;
          unsigned int xMax = (xBrick + 1) * xBrickDim_ - 1;
          unsigned int yMin = yBrick * yBrickDim_;
          unsigned int yMax = (yBrick + 1) * yBrickDim_ - 1;
          unsigned int zMin = zBrick * zBrickDim_;
          unsigned int zMax = (zBrick + 1) * zBrickDim_ - 1;
          unsigned int zLoc= 0;
          for (unsigned int zSub=zMin; zSub<=zMax; ++zSub) {
            unsigned int yLoc = 0;
            for (unsigned int ySub=yMin; ySub<=yMax; ++ySub) {  
              unsigned int xLoc = 0;
              for (unsigned int xSub=xMin; xSub<=xMax; ++xSub) {
                // Look up global index in full volume
                unsigned int globalIndex = 
                  xSub + ySub*xDim + zSub*xDim*yDim;
                // Set data at local subvolume index
                brick->SetData(xLoc, yLoc, zLoc, timestepData[globalIndex]); 
                xLoc++;
              }
              yLoc++;
            }
            zLoc++;
          }

          unsigned int brickIndex = 
            xBrick + yBrick*xNumBricks + zBrick*xNumBricks*yNumBricks;
          bricks_[brickIndex] = brick;
        
        }
      }
    }
  }

  return true;

}

bool Forge::Write() {
  return true;
}

