// This class defines some random number generation functions
// Most of the underlying code was borrowed from the BOOST library 
#ifndef __RANDOM_H__
#define __RANDOM_H__

#ifdef _MSC_VER
#include "stdintMSVC.h"
#else
#include <stdint.h>
#define uint64 uint64_t
#endif

#include <vector>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <limits>
#include <ctime>

namespace LaPSO {

namespace detail {

  template<bool is_signed>
  struct do_add
  { };

  template<>
  struct do_add<true>
  {
    template<class IntType>
    static IntType add(IntType m, IntType x, IntType c)
    {
      if (x < m - c)
        return x + c;
      else
        return x - (m-c);
    }
  };

  template<>
  struct do_add<false>
  {
    template<class IntType>
    static IntType add(IntType, IntType, IntType)
    {
      // difficult
      assert(!"const_mod::add with c too large");
      return 0;
    }
  };
//} // namespace detail


template<class IntType, IntType m>
class const_mod
{
public:
  static IntType add(IntType x, IntType c)
  {
    if(c == 0)
      return x;
    else if(c <= traits::max() - m)    // i.e. m+c < max
      return add_small(x, c);
    else
      return detail::do_add<traits::is_signed>::add(m, x, c);
  }

  static IntType mult(IntType a, IntType x)
  {
    if(a == 1)
      return x;
    else if(m <= traits::max()/a)      // i.e. a*m <= max
      return mult_small(a, x);
    else if(traits::is_signed && (m%a < m/a))
      return mult_schrage(a, x);
    else {
      // difficult
      assert(!"const_mod::mult with a too large");
      return 0;
    }
  }

  static IntType mult_add(IntType a, IntType x, IntType c)
  {
    if(m <= (traits::max()-c)/a)   // i.e. a*m+c <= max
      return (a*x+c) % m;
    else
      return add(mult(a, x), c);
  }

  static IntType invert(IntType x)
  { return x == 0 ? 0 : invert_euclidian(x); }

private:
  //typedef integer_traits<IntType> traits;
  typedef std::numeric_limits<IntType> traits;

  const_mod();      // don't instantiate

  static IntType add_small(IntType x, IntType c)
  {
    x += c;
    if(x >= m)
      x -= m;
    return x;
  }

  static IntType mult_small(IntType a, IntType x)
  {
    return a*x % m;
  }

  static IntType mult_schrage(IntType a, IntType value)
  {
    const IntType q = m / a;
    const IntType r = m % a;

    assert(r < q);        // check that overflow cannot happen

    value = a*(value%q) - r*(value/q);
    // An optimizer bug in the SGI MIPSpro 7.3.1.x compiler requires this
    // convoluted formulation of the loop (Synge Todo)
    for(;;) {
      if (value > 0)
        break;
      value += m;
    }
    return value;
  }

  // invert c in the finite field (mod m) (m must be prime)
  static IntType invert_euclidian(IntType c)
  {
    // we are interested in the gcd factor for c, because this is our inverse
    BOOST_STATIC_ASSERT(m > 0);
// #if BOOST_WORKAROUND(__MWERKS__, BOOST_TESTED_AT(0x3003))
//     assert(boost::integer_traits<IntType>::is_signed);
// #elif !defined(BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS)
//     BOOST_STATIC_ASSERT(boost::integer_traits<IntType>::is_signed);
// #endif
    assert(c > 0);
    IntType l1 = 0;
    IntType l2 = 1;
    IntType n = c;
    IntType p = m;
    for(;;) {
      IntType q = p / n;
      l1 -= q * l2;           // this requires a signed IntType!
      p -= q * n;
      if(p == 0)
        return (l2 < 1 ? l2 + m : l2);
      IntType q2 = n / p;
      l2 -= q2 * l1;
      n -= q2 * p;
      if(n == 0)
        return (l1 < 1 ? l1 + m : l1);
    }
  }
};

// The modulus is exactly the word size: rely on machine overflow handling.
// Due to a GCC bug, we cannot partially specialize in the presence of
// template value parameters.
template<>
class const_mod<unsigned int, 0>
{
  typedef unsigned int IntType;
public:
  static IntType add(IntType x, IntType c) { return x+c; }
  static IntType mult(IntType a, IntType x) { return a*x; }
  static IntType mult_add(IntType a, IntType x, IntType c) { return a*x+c; }

  // m is not prime, thus invert is not useful
private:                      // don't instantiate
  const_mod();
};

template<>
class const_mod<unsigned long, 0>
{
  typedef unsigned long IntType;
public:
  static IntType add(IntType x, IntType c) { return x+c; }
  static IntType mult(IntType a, IntType x) { return a*x; }
  static IntType mult_add(IntType a, IntType x, IntType c) { return a*x+c; }

  // m is not prime, thus invert is not useful
private:                      // don't instantiate
  const_mod();
};

// the modulus is some power of 2: rely partly on machine overflow handling
// we only specialize for rand48 at the moment
#ifndef _MSC_VER
template<>
class const_mod<uint64_t, uint64_t(1) << 48>
{
  typedef uint64_t IntType;
public:
  static IntType add(IntType x, IntType c) { return c == 0 ? x : mod(x+c); }
  static IntType mult(IntType a, IntType x) { return mod(a*x); }
  static IntType mult_add(IntType a, IntType x, IntType c)
    { return mod(a*x+c); }
  static IntType mod(IntType x) { return x &= ((uint64_t(1) << 48)-1); }

  // m is not prime, thus invert is not useful
private:                      // don't instantiate
  const_mod();
};
#endif
} // end namespace detail

// Because it is so commonly used: uniform distribution on the real [0..1)
// range.  This allows for specializations to avoid a costly int -> float
// conversion plus float multiplication
template<class UniformRandomNumberGenerator, class RealType = double>
class uniform_01
{
public:
  typedef UniformRandomNumberGenerator base_type;
  typedef RealType result_type;

  static const bool has_fixed_range = false;

//#if !defined(BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS) && !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300)
//  BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
//#endif

  explicit uniform_01(base_type rng)
    : _rng(rng),
      _factor(result_type(1) /
              (result_type((_rng.max)()-(_rng.min)()) +
               result_type(std::numeric_limits<base_result>::is_integer ? 1 : 0)))
  {
  }
  // compiler-generated copy ctor and copy assignment are fine

  result_type min() const { return result_type(0); }
  result_type max() const { return result_type(1); }
  base_type& base() { return _rng; }
  const base_type& base() const { return _rng; }
  void reset() { }

  result_type operator()() {
    for (;;) {
      result_type result = result_type(_rng() - (_rng.min)()) * _factor;
      if (result < result_type(1))
	return result;
    }
  }

//#if !defined(BOOST_NO_OPERATORS_IN_NAMESPACE) && !defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS)
//  template<class CharT, class Traits>
//  friend std::basic_ostream<CharT,Traits>&
//  operator<<(std::basic_ostream<CharT,Traits>& os, const uniform_01& u)
//  {
//    os << u._rng;
//    return os;
//  }
//
//  template<class CharT, class Traits>
//  friend std::basic_istream<CharT,Traits>&
//  operator>>(std::basic_istream<CharT,Traits>& is, uniform_01& u)
//  {
//    is >> u._rng;
//    return is;
//  }
//#endif

private:
  typedef typename base_type::result_type base_result;
  base_type _rng;
  result_type _factor;
};

// compile-time configurable linear congruential generator
template<class IntType, IntType a, IntType c, IntType m, IntType val>
class linear_congruential
{
private:    
  IntType _modulus;   // work-around for gcc "divide by zero" warning in ctor
  IntType _x;
public:
  typedef IntType result_type;
  static const bool has_fixed_range = true;
  static const result_type min_value = ( c == 0 ? 1 : 0 );
  static const result_type max_value = m-1;
  static const IntType multiplier = a;
  static const IntType increment = c;
  static const IntType modulus = m;

  explicit linear_congruential(IntType x0 = 1)
    : _modulus(modulus), _x(_modulus ? (x0 % _modulus) : x0)
  { 
    assert(c || x0); /* if c == 0 and x(0) == 0 then x(n) = 0 for all n */
    // overflow check
    // disabled because it gives spurious "divide by zero" gcc warnings
    // assert(m == 0 || (a*(m-1)+c) % m == (c < a ? c-a+m : c-a));
  }

  template<class It>
  linear_congruential(It& first, It last) { seed(first, last); }

  // compiler-generated copy constructor and assignment operator are fine
  void seed(IntType x0 = 1)
  {
    assert(c || x0);
    _x = (_modulus ? (x0 % _modulus) : x0);
  }

  template<class It>
  void seed(It& first, It last)
  {
    if(first == last)
      throw "linear_congruential::seed";
    IntType value = *first++;
    _x = (_modulus ? (value % _modulus) : value);
  }

  result_type min() const { return c == 0 ? 1 : 0; }
  result_type max() const { return modulus-1; }

  IntType operator()()
  {
    _x = detail::const_mod<IntType, m>::mult_add(a, _x, c);
    //_x = (a*_x+c) % m; // fails if a * _x > maximum integer
    return _x;
  }

  static bool validation(IntType x) { return val == x; }

  friend bool operator==(const linear_congruential& x,
                         const linear_congruential& y)
  { return x._x == y._x; }
  friend bool operator!=(const linear_congruential& x,
                         const linear_congruential& y)
  { return !(x == y); }
    
  template<class CharT, class Traits>
  friend std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT,Traits>& os,
             const linear_congruential& lcg)
  {
    return os << lcg._x;
  }

  template<class CharT, class Traits>
  friend std::basic_istream<CharT,Traits>&
  operator>>(std::basic_istream<CharT,Traits>& is,
             linear_congruential& lcg)
  {
    return is >> lcg._x;
  }
 
};
// validation values from the publications
typedef linear_congruential<int32_t, 16807, 0, 2147483647,1043618065> minstd_rand0;
typedef linear_congruential<int32_t, 48271, 0, 2147483647,399268537> minstd_rand;

  
template<class RealType, int w, unsigned int p, unsigned int q>
class lagged_fibonacci_01
{
public:
  typedef RealType result_type;
  static const bool  has_fixed_range = false; // or enum {has_fixed_range=false};
  static const int word_size = w;
  static const unsigned int long_lag = p;
  static const unsigned int short_lag = q;

  lagged_fibonacci_01() { init_modulus(); seed(); }
  explicit lagged_fibonacci_01(uint32_t value) { init_modulus(); seed(value); }
  template<class Generator>
  explicit lagged_fibonacci_01(Generator & gen) { init_modulus(); seed(gen); }
  template<class It> lagged_fibonacci_01(It& first, It last)
  { init_modulus(); seed(first, last); }
  // compiler-generated copy ctor and assignment operator are fine

private:
  void init_modulus()
  {
    _modulus = pow(RealType(2), word_size);
  }

public:
  void seed(uint32_t value = 331u)
  {
    minstd_rand0 intgen(value);
    seed(intgen);
  }

  // For GCC, moving this function out-of-line prevents inlining, which may
  // reduce overall object code size.  However, MSVC does not grok
  // out-of-line template member functions.
  template<class Generator>
  void seed(Generator & gen)
  {
     // use pass-by-reference, but wrap argument in pass_through_engine
     //typedef detail::pass_through_engine<Generator&> ref_gen;
     //uniform_01<ref_gen, RealType> gen01 =
     //  uniform_01<ref_gen, RealType>(ref_gen(gen));
	uniform_01<Generator, RealType> gen01(gen);
     // I could have used std::generate_n, but it takes "gen" by value
     for(unsigned int j = 0; j < long_lag; ++j)
       x[j] = gen01();
     i = long_lag;
  }

  template<class It>
  void seed(It& first, It last)
  {
    unsigned long mask = ~((~0u) << (w%32));   // now lowest w bits set
    RealType two32 = pow(RealType(2), 32);
    unsigned int j;
    for(j = 0; j < long_lag && first != last; ++j, ++first) {
      x[j] = RealType(0);
      for(int k = 0; k < w/32 && first != last; ++k, ++first)
        x[j] += *first / pow(two32,k+1);
      if(first != last && mask != 0)
        x[j] += fmod((*first & mask) / _modulus, RealType(1));
    }
    i = long_lag;
    if(first == last && j < long_lag)
      throw "lagged_fibonacci_01::seed";
  }

  result_type min() const { return result_type(0); }
  result_type max() const { return result_type(1); }

  result_type operator()()
  {
    if(i >= long_lag)
      fill();
    return x[i++];
  }

  static bool validation(result_type x)
  {
    throw "lagged_fibonacci_01::validation(x) not implemented";
    return false;
//     result_type v = fibonacci_validation<result_type, p, q>::value();
//     result_type epsilon = fibonacci_validation<result_type, p, q>::tolerance();
//     // std::abs is a source of trouble: sometimes, it's not overloaded
//     // for double, plus the usual namespace std noncompliance -> avoid it
//     // using std::abs;
//     // return abs(x - v) < 5 * epsilon
//     return x > v - epsilon && x < v + epsilon;
  }
  
//   template<class CharT, class Traits>
//   friend std::basic_ostream<CharT,Traits>&
//   operator<<(std::basic_ostream<CharT,Traits>& os, const lagged_fibonacci_01&f)
//   {
//     os << f.i << " ";
//     std::ios_base::fmtflags oldflags = os.flags(os.dec | os.fixed | os.left); 
//     for(unsigned int i = 0; i < f.long_lag; ++i)
//       os << f.x[i] * f._modulus << " ";
//     os.flags(oldflags);
//     return os;
//   }

//   template<class CharT, class Traits>
//   friend std::basic_istream<CharT, Traits>&
//   operator>>(std::basic_istream<CharT, Traits>& is, lagged_fibonacci_01& f)
//     {
//         is >> f.i >> std::ws;
//         for(unsigned int i = 0; i < f.long_lag; ++i) {
//             typename lagged_fibonacci_01::result_type value;
//             is >> value >> std::ws;
//             f.x[i] = value / f._modulus;
//         }
//         return is;
//     }

  friend bool operator==(const lagged_fibonacci_01& x,
                         const lagged_fibonacci_01& y)
  { return x.i == y.i && std::equal(x.x, x.x+long_lag, y.x); }
  friend bool operator!=(const lagged_fibonacci_01& x,
                         const lagged_fibonacci_01& y)
  { return !(x == y); }
#if 0
  // Use a member function; Streamable concept not supported.
  bool operator==(const lagged_fibonacci_01& rhs) const
  { return i == rhs.i && std::equal(x, x+long_lag, rhs.x); }
  bool operator!=(const lagged_fibonacci_01& rhs) const
  { return !(*this == rhs); }
#endif

private:
  void fill();
  unsigned int i;
  RealType x[long_lag];
  RealType _modulus;
};

template<class RealType, int w, unsigned int p, unsigned int q>
void lagged_fibonacci_01<RealType, w, p, q>::fill()
{
  // two loops to avoid costly modulo operations
  {  // extra scope for MSVC brokenness w.r.t. for scope
  for(unsigned int j = 0; j < short_lag; ++j) {
    RealType t = x[j] + x[j+(long_lag-short_lag)];
    if(t >= RealType(1))
      t -= RealType(1);
    x[j] = t;
  }
  }
  for(unsigned int j = short_lag; j < long_lag; ++j) {
    RealType t = x[j] + x[j-short_lag];
    if(t >= RealType(1))
      t -= RealType(1);
    x[j] = t;
  }
  i = 0;
}

typedef lagged_fibonacci_01<double, 48, 607, 273> lagged_fibonacci607;
typedef lagged_fibonacci_01<double, 48, 1279, 418> lagged_fibonacci1279;
typedef lagged_fibonacci_01<double, 48, 2281, 1252> lagged_fibonacci2281;
typedef lagged_fibonacci_01<double, 48, 3217, 576> lagged_fibonacci3217;
typedef lagged_fibonacci_01<double, 48, 4423, 2098> lagged_fibonacci4423;
typedef lagged_fibonacci_01<double, 48, 9689, 5502> lagged_fibonacci9689;
typedef lagged_fibonacci_01<double, 48, 19937, 9842> lagged_fibonacci19937;
typedef lagged_fibonacci_01<double, 48, 23209, 13470> lagged_fibonacci23209;
typedef lagged_fibonacci_01<double, 48, 44497, 21034> lagged_fibonacci44497;

class Uniform : public lagged_fibonacci1279 {
public:
  double operator()(double min=0.0,double max=1.0) {
    return min + (max-min)* ((*(lagged_fibonacci1279 *)this)());
  }
  void seedTime(){ // generate a random seed from current time
    seed((uint32_t)(time(0)%std::numeric_limits<uint32_t>::max()) );    
  }
};

}
#endif
