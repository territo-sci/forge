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

  std::cout << "Reading header" << std::endl;

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
  bricks_.resize(numTimesteps*xNumBricks*yNumBricks*zNumBricks);

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

  std::cout << "Reading volume data" << std::endl;

  // Read bricks
  for (unsigned int timestep=0; timestep<numTimesteps; ++timestep) {
  
    std::cout << "Reading timestep " << timestep << " of " << numTimesteps << 
              "\r" << std::flush;

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
          bricks_[timestep*xNumBricks*yNumBricks*zNumBricks+brickIndex]=brick;
        
        }
      }
    }
  }

  std::cout << "                                \r" << std::flush;
  std::cout << "Reading complete" << std::endl;

  return true;

}

bool Forge::Write() {

  std::fstream out;
  out.open(outFilename_.c_str(), std::ios_base::binary | std::ios_base::out |
           std::ios_base::trunc);
  if (out.fail()) {
    std::cout << "Failed to open " << outFilename_ << std::endl;
    return false;
  }

  unsigned int structure = header_->Structure();
  unsigned int dataDimensionality = header_->DataDimensionality();
  unsigned int xBrickDim = header_->XBrickDim();
  unsigned int yBrickDim = header_->YBrickDim();
  unsigned int zBrickDim = header_->ZBrickDim();
  unsigned int xNumBricks = header_->XNumBricks();
  unsigned int yNumBricks = header_->YNumBricks();
  unsigned int zNumBricks = header_->ZNumBricks();
  unsigned int numTimesteps = header_->NumTimesteps();
  unsigned int dataSize = header_->DataSize();

  out.seekp(std::ios_base::beg);

  std::cout << "Writing header" << std::endl;

  // Write header
  size_t s = sizeof(unsigned int);
  out.write(reinterpret_cast<char*>(&structure), s);
  out.write(reinterpret_cast<char*>(&dataDimensionality), s);
  out.write(reinterpret_cast<char*>(&xBrickDim), s);
  out.write(reinterpret_cast<char*>(&yBrickDim), s);
  out.write(reinterpret_cast<char*>(&zBrickDim), s);
  out.write(reinterpret_cast<char*>(&xNumBricks), s);
  out.write(reinterpret_cast<char*>(&yNumBricks), s);
  out.write(reinterpret_cast<char*>(&zNumBricks), s);
  out.write(reinterpret_cast<char*>(&numTimesteps), s);
  out.write(reinterpret_cast<char*>(&dataSize), s);

  std::cout << "Writing bricks" << std::endl;

  // Write bricks
  unsigned int i = 0;
  for (auto it=bricks_.begin(); it!=bricks_.end(); ++it) {
    std::cout << "Writing brick " << i++ << " of " << bricks_.size() << 
              "\r" << std::flush;
    out.write(reinterpret_cast<char*>(&((*it)->data_[0])), (*it)->Size());        
  }

  std::cout << "                                        \r" << std::flush;
  std::cout << "Writing complete" << std::endl;

  out.close();

  return true;
}

