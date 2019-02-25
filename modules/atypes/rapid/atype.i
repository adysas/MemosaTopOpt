%{
#include "Rapid.h"
  %}

#if 0
%typemap(in) Rapid
{
  if (PyFloat_Check($input))
  {
      $1 = Rapid(PyFloat_AsDouble($input));
  }
  else if (PyInt_Check($input))
  {
      $1 = Rapid(double(PyInt_AsLong($input)));
  }
  else if (PyTuple_Check($input) && (2==PyTuple_Size($input)))
  {
      $1 = Rapid(PyFloat_AsDouble(PyTuple_GetItem($input,0)),
                   PyFloat_AsDouble(PyTuple_GetItem($input,1)));
  }
  else
    throw CException("invalid Rapid input");
}
#endif

class Rapid
{
public:
  typedef map<int, double> PartialDerivatives;
  Rapid();
  Rapid(double a, PartialDerivatives& b);
  double _v;
  PartialDerivatives _dv;
};

typedef Rapid ATYPE;


//%include "atype.swg"

#define USING_ATYPE_RAPID 

#define ATYPE_STR Rapid

