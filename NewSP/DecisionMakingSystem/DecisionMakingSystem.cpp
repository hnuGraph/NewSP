#include "DecisionMakingSystem.h"
//adaptive check
LRAndIndexCheckType DMS::makeDecision(const vertexType & currentType, const uint & destListSize, const uint & unfreezeListSize){
    if(currentType == isolatedVertex){//isolated not check
        return ioslatedVertexNothing;
    }
    if(unfreezeListSize == 0){//do not need to unfreeze
        if(currentType == freeVertex){
            return Part1Check;//just freeVertex need to check
        }
        return Part1Nothing;
    }
    else{/*need to unfreeze*/
        if(destListSize != 0){
            return Part1Check;//if have expend dependent, just check
        }
        else{
            if(currentType == freeVertex){
                return Part2Check;//just freeVertex need to check
            }
            return Part2Nothing;
        }
    }
}

