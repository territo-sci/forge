/*
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#include <BricksHeader.h>

using namespace osp;

BricksHeader * BricksHeader::New() {
  return new BricksHeader();
}

BricksHeader::BricksHeader()
  : structure_(0), dataDimensionality_(0), xBrickDim_(0), yBrickDim_(0),
    zBrickDim_(0), xNumBricks_(0), yNumBricks_(0), zNumBricks_(0),
    numTimesteps_(0), dataSize_(0) {
}

BricksHeader::~BricksHeader() {
}

void BricksHeader::SetStructure(unsigned int _structure) {
  structure_ = _structure;
}

void BricksHeader::SetDataDimensionality(unsigned int _dataDimensionality) {
  dataDimensionality_ = _dataDimensionality;
}

void BricksHeader::SetDimensions(unsigned int _xBrickDim,
                                 unsigned int _yBrickDim,
                                 unsigned int _zBrickDim) {
  xBrickDim_ = _xBrickDim;
  yBrickDim_ = _yBrickDim;
  zBrickDim_ = _zBrickDim;
}

void BricksHeader::SetNumBricks(unsigned int _xNumBricks,
                                unsigned int _yNumBricks,
                                unsigned int _zNumBricks) {
  xNumBricks_ = _xNumBricks;
  yNumBricks_ = _yNumBricks;
  zNumBricks_ = _zNumBricks;
}

void BricksHeader::SetNumTimesteps(unsigned int _numTimesteps) {
  numTimesteps_ = _numTimesteps;
}

void BricksHeader::SetDataSize(unsigned int _dataSize) {
  dataSize_ = _dataSize;
}


