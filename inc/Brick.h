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
  static Brick * New(size_t _xDim,
                     size_t _yDim,
                     size_t _zDim,
                     T _defaultVal) {
    return new Brick<T>(_xDim, _yDim, _zDim, _defaultVal);
  }
  ~Brick() {}
  
  void SetData(size_t _x,
               size_t _y,
               size_t _z,
               T _value) {
    data_[Index(_x, _y, _z)] = _value;
  }

  T Data(size_t _x,
         size_t _y,
         size_t _z) const {
    return data_[Index(_x, _y, _z)];
  }

  size_t Size() const {
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
    for (size_t z=0; z<zDim_; ++z) {
      for (size_t y=0; y<yDim_; ++y) {
        for (size_t x=0; x<xDim_; ++x) {
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
    for (size_t z=0; z<zDim_; ++z) {
      for (size_t y=0; y<yDim_; ++y) {
        for (size_t x=0; x<xDim_; ++x) {
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

  // Filter a brick by combining 2x2x2 voxels
  // Resulting brick has half size dimensions
  static Brick<T> * Filter(Brick<T> *_brick) {
    size_t dim = _brick->xDim_;
    Brick<T> *out = Brick<T>::New(dim/2, dim/2, dim/2, static_cast<T>(0));
    for (size_t z=0; z<dim; z+=2) {
      for (size_t y=0; y<dim; y+=2) {
        for (size_t x=0; x<dim; x+=2) {
          std::vector<T> toAverage(8, static_cast<T>(0));
          // Specify the eight voxels to average
          toAverage[0] = _brick->Data(x,   y,   z  );
          toAverage[1] = _brick->Data(x+1, y,   z  );
          toAverage[2] = _brick->Data(x,   y+1, z  );
          toAverage[3] = _brick->Data(x+1, y+1, z  );
          toAverage[4] = _brick->Data(x,   y,   z+1);
          toAverage[5] = _brick->Data(x+1, y,   z+1);
          toAverage[6] = _brick->Data(x,   y+1, z+1);
          toAverage[7] = _brick->Data(x+1, y+1, z+1);
          T sum = static_cast<T>(0);
          for (auto it=toAverage.begin(); it!=toAverage.end(); ++it) {
            sum += *it;
          }
          out->SetData(x/2, y/2, z/2, sum/static_cast<T>(8));
        }
      }
    }
    return out;
  }

  // Combine eight bricks and return a new one, where the dimensions
  // are twice the dimensions of the combined bricks.
  // The combined bricks should be in a vector, ordered as follows:
  // 0: x,   y,   z
  // 1: x+1, y,   z
  // 2: x,   y+1, z
  // 3: x+1, y+1, z
  // 4: x,   y,   z+1
  // 5: x+1, y,   z+1
  // 6: x,   y+1, z+1
  // 7: x+1, y+1, z+1
  static Brick<T> * Combine(std::vector<Brick<T> *> _bricks) {
    // Assume all input sizes are equal
    size_t dim = _bricks[0]->xDim_;
    Brick<T> *out =  Brick<T>::New(dim*2, dim*2, dim*2, static_cast<T>(0));
    
    if (_bricks.size() != 8) {
      std::cout << "ERROR: Combine vector not of size 8" << std::endl;
      return out;
    }

    // Loop over positions in new brick 
    for (size_t z=0; z<dim*2; ++z) {
      for (size_t y=0; y<dim*2; ++y) {
        for (size_t x=0; x<dim*2; ++x) {
        
          T value;
          // Choose and sample from the right quadrant
          if (x < dim) { // quadrant 0, 2, 4 or 6
            if (y < dim) { // quadrant 0 or 4
              if (z < dim) { // quadrant 0
                value = _bricks[0]->Data(x, y, z); 
              } else { // quadrant 4
                value = _bricks[4]->Data(x, y, z-dim);
              }
            } else { // quadrant 2 or 6
              if (z < dim) { // quadrant 2
                value = _bricks[2]->Data(x, y-dim, z);
              } else { // quadrant 6
                value = _bricks[6]->Data(x, y-dim, z-dim);
              }
            }
          } else {// quadrant 1, 3, 5 or 7
            if (y < dim) { // quadrant 1 or 5
              if (z < dim) { // quadrant 1
                value = _bricks[1]->Data(x-dim, y, z);
              } else { // quadrant 5
                value = _bricks[5]->Data(x-dim, y, z-dim); 
              } 
            } else { // quadrant 3 or 7
              if (z < dim) { // quadrant 3
                value = _bricks[3]->Data(x-dim, y-dim, z);
              } else { // quadrant 7
                value = _bricks[7]->Data(x-dim, y-dim, z-dim);
              }
            }
          }

          out->SetData(x, y, z, value);

        }
      }
    }

    return out;

  }

  friend class Forge;

private:
  Brick(size_t _xDim,
        size_t _yDim,
        size_t _zDim,
        T _defaultVal)
        : xDim_(_xDim), yDim_(_yDim), zDim_(_zDim) {
    data_.resize(_xDim*_yDim*_zDim, _defaultVal);
  }
  Brick();
  Brick(const Brick&);
  size_t Index(size_t _x,
                     size_t _y,
                     size_t _z) const {
    return _x + _y*xDim_ + _z*xDim_*yDim_;
  }
  size_t xDim_;
  size_t yDim_;
  size_t zDim_;
  std::vector<T> data_;
};

}


