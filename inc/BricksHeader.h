/*
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#ifndef BRICKSHEADER_H_
#define BRICKSHEADER_H_

namespace osp {

class BricksHeader {
public:
  static BricksHeader * New();
  ~BricksHeader();
  
  unsigned int Structure() const { return structure_; }
  unsigned int DataDimensionality() const { return dataDimensionality_; }
  unsigned int XBrickDim() const { return xBrickDim_; }
  unsigned int YBrickDim() const { return yBrickDim_; }
  unsigned int ZBrickDim() const { return zBrickDim_; }
  unsigned int XNumBricks() const { return xNumBricks_; }
  unsigned int YNumBricks() const { return yNumBricks_; }
  unsigned int ZNumBricks() const { return zNumBricks_; }
  unsigned int NumTimesteps() const { return numTimesteps_; }
  unsigned int DataSize() const { return dataSize_; }
    
  void SetStructure(unsigned int _structure);
  void SetDataDimensionality(unsigned int _dataDimensionality);
  void SetDimensions(unsigned int _xDim,
                     unsigned int _yDim,
                     unsigned int _zDim);
  void SetNumBricks(unsigned int _xNumBricks,
                    unsigned int _yNumBricks,
                    unsigned int _zNumBricks);
  void SetNumTimesteps(unsigned int _numTimesteps);
  void SetDataSize(unsigned int _dataSize);

private:
  BricksHeader();
  BricksHeader(const BricksHeader&);

  // Encodes the structure
  // 0 - flat structure
  // TODO 1 - octree
  // TODO 2 - tsp?
  unsigned int structure_;
  // Data dimensionality (1 for scalar, 3 for 3D vector etc)
  unsigned int dataDimensionality_;
  // Dimensions for individual bricks
  unsigned int xBrickDim_;
  unsigned int yBrickDim_;
  unsigned int zBrickDim_;
  // Number of bricks per axis (for flat structure)
  unsigned int xNumBricks_;
  unsigned int yNumBricks_;
  unsigned int zNumBricks_;
  // Number of timesteps in file
  unsigned int numTimesteps_;
  // Size of brick data entry, in bytes
  unsigned int dataSize_;
  
};

}

#endif

