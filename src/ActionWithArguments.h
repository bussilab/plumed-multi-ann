#ifndef __PLUMED_ActionWithArguments_h
#define __PLUMED_ActionWithArguments_h

#include "Action.h"
#include <map>
#include <vector>
#include <cassert>

#include "ActionWithValue.h"

namespace PLMD{

/// Action which takes other Action's as arguments.
/// Arguments are objects of type PLMD::Value, and
/// are addressed using the ARG= keyword on the directive line
class ActionWithArguments:
  public virtual Action
{
  std::vector<Value*> arguments;

protected:
                           ActionWithArguments(const ActionOptions&);
  virtual                 ~ActionWithArguments(){};
  std::vector<Value*>    & getArguments();
  double                   getArgument(int)const;
  unsigned                 getNumberOfArguments()const;
public:
};


inline
std::vector<Value*> & ActionWithArguments::getArguments(){
  return arguments;
}

inline
double ActionWithArguments::getArgument(int i)const{
  return arguments[i]->getValue();
}

inline
unsigned ActionWithArguments::getNumberOfArguments()const{
  return arguments.size();
}


}

#endif
