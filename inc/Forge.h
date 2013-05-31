/*
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#ifndef FORGE_H_
#define FORGE_H_

#include <string>
#include <vector>

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
  void SetStructure(unsigned int structure_);
  void SetBrickDimensions(unsigned int _xBrickDim,
                          unsigned int _yBrickDim,
                          unsigned int _zBrickDim);

  bool Read();
  bool Write();

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

};


}

#endif
