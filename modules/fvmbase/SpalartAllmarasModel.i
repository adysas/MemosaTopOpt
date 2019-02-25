%{
#include "SpalartAllmarasModel.h"
%}

using namespace std;

%include "FloatVarDict.i"
%include "Model.i"

template <class T>
struct SpalartAllmarasBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct SpalartAllmarasVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct SpalartAllmarasModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  bool useCentralDifference;
  bool transient;
}; 


%template(SpalartAllmarasBCA) SpalartAllmarasBC< ATYPE_STR >;
%template(SpalartAllmarasVCA) SpalartAllmarasVC< ATYPE_STR >;
%template(SpalartAllmarasBCList) std::vector<SpalartAllmarasBC< ATYPE_STR >* >;
%template(SpalartAllmarasBCsMap) std::map<int,SpalartAllmarasBC< ATYPE_STR >* >;
%template(SpalartAllmarasVCList) std::vector<SpalartAllmarasVC< ATYPE_STR >* >;
%template(SpalartAllmarasVCsMap) std::map<int,SpalartAllmarasVC< ATYPE_STR >* >;

%template(SpalartAllmarasModelOptionsA) SpalartAllmarasModelOptions< ATYPE_STR >;


%import "Model.i"

%include "SpalartAllmarasModel.h"


%template(SpalartAllmarasModelA) SpalartAllmarasModel< ATYPE_STR >;

