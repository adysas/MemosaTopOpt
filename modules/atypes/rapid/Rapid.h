#ifndef _RAPID_H_
#define _RAPID_H_
#include "NumType.h"
#include "misc.h"

#include <string.h>
#include <map>
#include <utility>
#include <boost/foreach.hpp>
#include "math.h"
using namespace std;
using namespace boost;
#define foreach         BOOST_FOREACH


class Rapid
{
public:
  typedef map<int, double> PartialDerivatives;
  typedef map<int, double>::iterator Iter;
  typedef map<int, double>::const_iterator CIter;

  typedef Rapid This_T;
  typedef Rapid T_Scalar;
  typedef NumTypeTraits<double>::T_BuiltIn T_BuiltIn;

  static string getTypeName()
  {
    return "Rapid";
  }
  
  // no initialization for empty ctor - consistent with builtins

  ///////////////////////////////////////////////////////////////////
  Rapid()
  {}
  
  Rapid(const double v, const int key):
    _v(v)
  {
    _dv[key] = 1;
  }
  
  Rapid(const double v, const PartialDerivatives& dv) :
    _v(v),
    _dv(dv)
  {}
  
  //Used mostly for setting independent parameters
  explicit Rapid(const double v) :
    _v(v)
  {}
  
  explicit Rapid(const int v) :
    _v(double(v))
  {}

  //Copy constructor
  Rapid(const Rapid& o) :
    _v(o._v),
    _dv(o._dv)
  {}

  ~Rapid()
  {}

  //////////////////////////////////////////////////////////////////////
  void setIndex(const int& key)
  {
    _dv[key] = 1;
  }
    
  //////////////////////////////////////////////////////////////////////
  Rapid& operator=(const Rapid& o)
  {
    if (this == &o)
      return *this;
    _v = o._v;
    _dv = o._dv;
    return *this;
  }

  Rapid& operator=(const double& f)
  {
    _v = f;
    return *this;
  }

  Rapid& operator=(const int& i)
  {
    _v = double(i);     
    return *this;
  }
  //////////////////////////////////////////////////////////////////////
  Rapid& operator+=(const Rapid& o)
  {
    _v += o._v;
    
    foreach(const typename PartialDerivatives::value_type& ii, o._dv)
      {
	_dv[ii.first] += ii.second;
      }
    
    return *this;
  }

  Rapid& operator-=(const Rapid& o)
  {
    _v -= o._v;
    
    foreach(const typename PartialDerivatives::value_type& ii, o._dv)
      {
	_dv[ii.first] -= ii.second;
      }
    
    return *this;
  }

  Rapid& operator*=(const Rapid& o)
  {

    foreach(const typename PartialDerivatives::value_type& ii, o._dv)
      {
	Iter key;
	key = _dv.find(ii.first);
	if (key != _dv.end())	
	  {
	    _dv[ii.first] = _dv[ii.first]*o._v + _v*ii.second;
	  }
	else
	  {
	    _dv[ii.first] = _v*ii.second;
	  }
      }
    
    foreach(const typename PartialDerivatives::value_type& ii, _dv)
      {
	CIter key;
	key = o._dv.find(ii.first);
	if (key == o._dv.end())	
	  {
	    _dv[ii.first] = _dv[ii.first]*o._v;
	  }
      }
          
    _v *= o._v;
    return *this;
  }


  Rapid& operator/=(const Rapid& o)
  {

    foreach(const typename PartialDerivatives::value_type& ii, o._dv)
      {
	Iter key;
	key = _dv.find(ii.first);
	if (key != _dv.end())	
	  {
	    _dv[ii.first] = (_dv[ii.first]*o._v - _v*ii.second)/(o._v*o._v);
	  }
	else
	  {
	    _dv[ii.first] = (-1*_v*ii.second)/(o._v*o._v);
	  }
      }
    
    foreach(const typename PartialDerivatives::value_type& ii, _dv)
      {
	CIter key;
	key = o._dv.find(ii.first);
	if (key == o._dv.end())	
	  {
	    _dv[ii.first] = (_dv[ii.first]*o._v)/(o._v*o._v);
	  }
      }
          
    _v /= o._v;
    return *this;
  }


  //////////////////////////////////////////////////////////////////////
  Rapid& operator+=(const double& o)
  {
    _v += o;
    return *this;
  }

  Rapid& operator-=(const double& o)
  {
    _v -= o;
    return *this;
  }

  Rapid& operator*=(const double& o)
  {
    _v *= o;
      
    foreach(const typename PartialDerivatives::value_type& ii, _dv)
      {
	_dv[ii.first] *= o;
      }
    return *this;
  }

  Rapid& operator/=(const double& o)
  {
    _v /= o;
      
    foreach(const typename PartialDerivatives::value_type& ii, _dv)
      {
	_dv[ii.first] /= o;
      }
    return *this;
  }

  //////////////////////////////////////////////////////////////////////
  Rapid& operator+=(const int& i)
  {
    _v += double(i);
    return *this;
  }

  Rapid& operator-=(const int& i)
  {
    _v -= double(i);
    return *this;
  }

  Rapid& operator*=(const int& i)
  {
    _v *= double(i);
      
    foreach(const typename PartialDerivatives::value_type& ii, _dv)
      {
	_dv[ii.first] *= double(i);	
      }
    return *this;
  }

  Rapid& operator/=(const int& i)
  {
    _v /= double(i);
      
    foreach(const typename PartialDerivatives::value_type& ii, _dv)
      {
	_dv[ii.first] /= double(i);
      }
    return *this;
  }
  //////////////////////////////////////////////////////////////////////


#define RAPID_RELATIONAL_OPS(opname,_op_)	\
  bool opname(const Rapid& o) const		\
  {                                             \
    return (_v _op_ o._v);                      \
  }                                             \
  bool opname(const double& o) const		\
  {                                             \
    return (_v _op_ o);                         \
  }                                             \
  bool opname(const int& o) const		\
  {                                             \
    return (_v _op_ double(o));                 \
  }

  RAPID_RELATIONAL_OPS(operator>,>);
  RAPID_RELATIONAL_OPS(operator>=,>=);
  RAPID_RELATIONAL_OPS(operator<,<);
  RAPID_RELATIONAL_OPS(operator<=,<=);
  RAPID_RELATIONAL_OPS(operator==,==);
  RAPID_RELATIONAL_OPS(operator!=,!=);
  
#undef RAPID_RELATIONAL_OPS

  /////////////////////////////////////////////////////////////////////

  /*
  void print(std::ostream &os) const
  {
    os << "< "  << _v << ": "; 
    
    foreach(const typename PartialDerivatives::value_type& ii, _dv)
      {
	os << ii.first << ", " << ii.second << "; ";
      }
    os << "> "; 
  }
  */
  void print(std::ostream& os) const 
  {
    os << "<val = " << _v;   
    foreach(const  typename PartialDerivatives::value_type& ii, _dv)
      {
	os << ", \u2202[" << ii.first << "] = " << ii.second;
      }
    os << "> ";
  }


  static Rapid getZero()
  {
    double zero = NumTypeTraits<double>::getZero();
    return Rapid(zero);
  }

  static Rapid getUnity()
  {
    return Rapid(NumTypeTraits<double>::getUnity());
  }

  static Rapid getNegativeUnity()
  {
    return Rapid(NumTypeTraits<double>::getNegativeUnity());
  }

  static double doubleMeasure(const Rapid& x)
  {
    return NumTypeTraits<double>::doubleMeasure(x._v);
  }
  
  /*
  static void setFloat(Rapid& t, const int i, const double& val)
  {
    throw;
    if (i == 0)
      t._v = val;
    else
    t._dv = val;
  }

  static double getFloat(const Rapid& t, const int i)
  {
    throw;
    if (i==0)
      return t._v;
    else
    return t._dv;
  }*/
  
  // only printing the value and not the derivative
  //static void write(FILE* fp, const Tangent& x) {fprintf(fp,"%f",x._v);}

  static int getDimension() 
  {
    //return NumTypeTraits<double>::getDimension()+1;
    throw;
    return 0;
  }
  
  static void getShape(int *shp) 
  { 
    //*shp = 2; 
    //NumTypeTraits<double>::getShape(shp+1);
    throw;
  }

  
  static int getDataSize() 
  {
    throw;
    //return 0;
  }
  

  Rapid fabs() const
  {
    Rapid o;
    o._v = ::fabs(_v);
    
    foreach(const typename Rapid::PartialDerivatives::value_type& ii, _dv)
    {
      o._dv[ii.first] = (_v > 0.0 ? ii.second : -ii.second); 
    }
  return o;
  }
    
  static void accumulateOneNorm(Rapid& sum, const Rapid& v) 
  { 
    throw;
  }
  
  static void accumulateDotProduct(Rapid& sum, const Rapid& v0, const Rapid& v1)
  {
    throw;
  }
  
  static void reduceSum(Rapid& sum, const Rapid& x) {sum+=x;}

  static void safeDivide(Rapid& x, const Rapid& y) {if (y._v!=0) x/=y;}
  static void normalize(Rapid& x, const Rapid& y) {if (y._v!=0) x/=y;}
  static void setMax(Rapid& x, const Rapid& y) {if (y._v>x._v) x=y;}
  
  
  void clearRapid()
  {
    _dv.clear();
  }
    
  double _v;
  PartialDerivatives _dv;
};



#define RAPID_BINARY_OP(opname,_op_)		\
  Rapid opname(const Rapid& a, const Rapid& b)	\
  {						\
    return Rapid(a) _op_ b;			\
  }						\
  Rapid opname(const Rapid& a, const double& b)	\
  {						\
    return Rapid(a) _op_ b;			\
  }						\
  Rapid opname(const double& a, const Rapid& b)	\
  {						\
    return Rapid(a) _op_ b;			\
  }						\
  Rapid opname(const Rapid& a, const int& b)	\
  {						\
    return Rapid(a) _op_ b;			\
  }						\
  Rapid opname(const int& a, const Rapid& b)	\
  {						\
    return Rapid(a) _op_ b;			\
  }						\


RAPID_BINARY_OP(operator+,+=);
RAPID_BINARY_OP(operator*,*=);
RAPID_BINARY_OP(operator-,-=);
RAPID_BINARY_OP(operator/,/=);



Rapid operator+(const Rapid& a)
{
  return a;
}

Rapid operator-(const Rapid& a)
{
  Rapid o;
  o._v = -a._v; 
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = -ii.second;      
    }
  return o;
}

inline std::ostream &operator<<(std::ostream &os, const Rapid &a)
{
  a.print(os);
  return os;
}

Rapid cos(const Rapid& a)
{
  Rapid o;
  o._v = cos(a._v);
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = ii.second*-sin(ii.first);
    }
  return o;
}

Rapid sin(const Rapid& a)
{
  Rapid o;
  o._v = sin(a._v);
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = ii.second*cos(ii.first);
    }
  return o;
}

Rapid exp(const Rapid& a)
{
  Rapid o;
  o._v = exp(a._v);
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = ii.second*exp(ii.first);
    }
  return o;
}


Rapid pow(const Rapid& a, const double& b)
{
  Rapid o;
  o._v = a._v == 0.0 ? 0 : pow(a._v,b); 
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = a._v == 0.0 ? 0 : ii.second*b*pow(a._v,b-1);      
    }
  return o;
}

Rapid pow(const Rapid& a, const Rapid& b)
{
  Rapid o;
  o._v = a._v == 0.0 ? 0 : pow(a._v,b._v); 
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = a._v == 0.0 ? 0 : ii.second*b._v*pow(a._v,b._v-1);      
    }
  return o;
}

Rapid abs(const Rapid& a)
{
  Rapid o;
  o._v = ::fabs(a._v);
    
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = (a._v > 0.0 ? ii.second : -ii.second); 
    }
  return o;
}

Rapid fabs(const Rapid& a)
{
  Rapid o;
  o._v = fabs(a._v);
    
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = (a._v > 0.0 ? ii.second : -ii.second); 
    }
  return o;
}

Rapid sqrt(const Rapid& a)
{
  Rapid o;
  double sqv = sqrt(a._v);
  o._v = sqv;
  foreach(const typename Rapid::PartialDerivatives::value_type& ii, a._dv)
    {
      o._dv[ii.first] = (sqv == 0.0 ? 0 : 0.5*ii.second/sqv); 
    }
  return o;
}

double ceil(const Rapid& a)
{
  return ceil(a._v);
}


#endif
