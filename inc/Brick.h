/* 
 * Author: Victor Sand (victor.sand@gmail.com)
 *
 */

#include <iostream>
#include <vector>

namespace osp {

template <class T>
class Brick {
public:
  static Brick * New(unsigned int _xDim,
                     unsigned int _yDim,
                     unsigned int _zDim,
                     T _defaultVal) {
    return new Brick<T>(_xDim, _yDim, _zDim, _defaultVal);
  }
  ~Brick() {}
  
  void SetData(unsigned int _x,
               unsigned int _y,
               unsigned int _z,
               T _value) {
    data_[Index(_x, _y, _z)] = _value;
  }

  T Data(unsigned int _x,
         unsigned int _y,
         unsigned int _z) const {
    return data_[Index(_x, _y, _z)];
  }

  unsigned int Size() const {
    return sizeof(T)*data_.size();
  }

  friend class Forge;

private:
  Brick(unsigned int _xDim,
        unsigned int _yDim,
        unsigned int _zDim,
        T _defaultVal)
        : xDim_(_xDim), yDim_(_yDim), zDim_(_zDim) {
    data_.resize(_xDim*_yDim*_zDim, _defaultVal);
  }
  Brick();
  Brick(const Brick&);
  unsigned int Index(unsigned int _x,
                     unsigned int _y,
                     unsigned int _z) const {
    return _x + _y*xDim_ + _z*xDim_*yDim_;
  }
  unsigned int xDim_;
  unsigned int yDim_;
  unsigned int zDim_;
  std::vector<T> data_;
};

}


