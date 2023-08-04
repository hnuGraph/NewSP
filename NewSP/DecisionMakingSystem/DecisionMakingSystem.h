#ifndef DecisionMakingSystem_DecisionMakingSystem
#define DecisionMakingSystem_DecisionMakingSystem
#include "utils/types.h"
#include "utils/globals.h"
class DMS{
public:
    LRAndIndexCheckType makeDecision(const vertexType & currentType, const uint & destListSize, const uint & unfreezeListSize);
};
#endif