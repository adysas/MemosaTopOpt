%{
#include "WallDistanceModel.h"
%}

using namespace std;

%include "FloatVarDict.i"
%include "Model.i"

template <class T>
struct WallDistanceBC : public FloatVarDict<T>
{
  string bcType;
}; 


template <class T>
struct WallDistanceModelOptions : public FloatVarDict<T>
{
  double relativeTolerance;
  double absoluteTolerance;
  LinearSolver *linearSolver;
}; 


%template(WallDistanceBCA) WallDistanceBC< ATYPE_STR >;
%template(WallDistanceBCList) std::vector<WallDistanceBC< ATYPE_STR >* >;
%template(WallDistanceBCsMap) std::map<int,WallDistanceBC< ATYPE_STR >* >;

%template(WallDistanceModelOptionsA) WallDistanceModelOptions< ATYPE_STR >;


%import "Model.i"

%include "WallDistanceModel.h"


%template(WallDistanceModelA) WallDistanceModel< ATYPE_STR >;

