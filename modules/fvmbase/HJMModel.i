%{
#include "HamiltonJacobiModel.h"
%}

using namespace std;

%include "FloatVarDict.i"
%include "Model.i"

template <class T>
struct HamiltonJacobiBC : public FloatVarDict<T>
{
  string bcType;
}; 

template <class T>
struct HamiltonJacobiVC : public FloatVarDict<T>
{
  string vcType;
}; 

template <class T>
struct HamiltonJacobiModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
  bool useCentralDifference;
}; 


%template(HamiltonJacobiBCA) HamiltonJacobiBC< ATYPE_STR >;
%template(HamiltonJacobiVCA) HamiltonJacobiVC< ATYPE_STR >;
%template(HamiltonJacobiBCList) std::vector<HamiltonJacobiBC< ATYPE_STR >* >;
%template(HamiltonJacobiBCsMap) std::map<int,HamiltonJacobiBC< ATYPE_STR >* >;
%template(HamiltonJacobiVCList) std::vector<HamiltonJacobiVC< ATYPE_STR >* >;
%template(HamiltonJacobiVCsMap) std::map<int,HamiltonJacobiVC< ATYPE_STR >* >;

%template(HamiltonJacobiModelOptionsA) HamiltonJacobiModelOptions< ATYPE_STR >;


%import "Model.i"

%include "HamiltonJacobiModel.h"


%template(HamiltonJacobiModelA) HamiltonJacobiModel< ATYPE_STR >;

