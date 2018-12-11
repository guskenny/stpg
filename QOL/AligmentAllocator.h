#ifndef ALIGNMENT_ALLOCATOR_H
#define ALIGNMENT_ALLOCATOR_H

//if needed, a vector can be declared like this
//std::vector<double, AlignmentAllocator<double, 32> > dd;

#include <stdlib.h>
#include <malloc.h>

template <typename T, size_t N = 16>
class AlignmentAllocator {
public:
  typedef T value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  typedef T * pointer;
  typedef const T * const_pointer;

  typedef T & reference;
  typedef const T & const_reference;

  public:
  inline AlignmentAllocator () throw () { }

  template <typename T2>
  inline AlignmentAllocator (const AlignmentAllocator<T2, N> &) throw () { }

  inline ~AlignmentAllocator () throw () { }

  inline pointer adress (reference r) {
    return &r;
  }

  inline const_pointer adress (const_reference r) const {
    return &r;
  }

  inline pointer allocate (size_type n) {
#ifdef _MSC_VER
     return (pointer)_aligned_malloc(n*sizeof(value_type), N);
#else
    return (pointer)__mm_malloc(n*sizeof(value_type), N);
#endif
  }

  inline void deallocate (pointer p, size_type) {
#ifdef _MSC_VER
    _aligned_free (p);
#else
    mm_free(p);
#endif
  }

  inline void construct (pointer p, const value_type & wert) {
     new (p) value_type (wert);
  }

  inline void destroy (pointer p) {
    p->~value_type ();
  }

  inline size_type max_size () const throw () {
    return size_type (-1) / sizeof (value_type);
  }

  template <typename T2>
  struct rebind {
    typedef AlignmentAllocator<T2, N> other;
  };

  bool operator!=(const AlignmentAllocator<T,N>& other) const  {
    return !(*this == other);
  }

  // Returns true if and only if storage allocated from *this
  // can be deallocated from other, and vice versa.
  // Always returns true for stateless allocators.
  bool operator==(const AlignmentAllocator<T,N>& other) const {
    return true;
  }
};

#endif