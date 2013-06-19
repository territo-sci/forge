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

  // Add the values of one brick
  bool Add(Brick<T> *_brick) {
    if (xDim_ != _brick->xDim_ ||
        yDim_ != _brick->yDim_ ||
        zDim_ != _brick->zDim_) {
      std::cout << "Brick::Add(): Dimension mismatch" << std::endl;
      return false;
    }
    for (unsigned int z=0; z<zDim_; ++z) {
      for (unsigned int y=0; y<yDim_; ++y) {
        for (unsigned int x=0; x<xDim_; ++x) {
          SetData(x, y, z, Data(x, y, z) + _brick->Data(x, y, z));  
        }
      }
    }
    return true;
  }

  // Divide all brick values by a number
  bool Divide(T _divisor) {
    if (_divisor == 0) {
      std::cout << "Brick::Divide(): Division by zero" << std::endl;
      return false;
    }
    for (unsigned int z=0; z<zDim_; ++z) {
      for (unsigned int y=0; y<yDim_; ++y) {
        for (unsigned int x=0; x<xDim_; ++x) {
          SetData(x, y, z, Data(x, y, z)/_divisor);  
        }
      }
    }
    return true;
  }

  // Average two bricks, return result
  static Brick<T> * Average(Brick<T> *_a, Brick<T> *_b) {
    if (_a->xDim_ != _b->xDim_ ||
        _a->yDim_ != _b->yDim_ ||
        _a->zDim_ != _b->zDim_) {
      std::cout << "Brick::Average(): Dimension mismatch" << std::endl;
      return NULL;
    }
    Brick<T> *out = New(_a->xDim_, _a->yDim_, _a->zDim_, static_cast<T>(0));
    out->Add(_a);
    out->Add(_b);
    out->Divide(static_cast<T>(2));
    return out;
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


