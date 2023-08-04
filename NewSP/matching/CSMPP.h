#include <vector>
#include <map>

#include "utils/types.h"
#include "utils/globals.h"
#include "graph/graph.h"
#include "matching/matching.h"

enum searchType{
    init, pos, neg
};
class CSMPP : public matching
{
private:
    
    std::vector<std::pair<Edge, std::vector<std::pair<uint,uint>>>> updateEdgeFindQueryEqual;// record the edge update 
    std::map<Edge,std::vector<std::pair<uint,uint>>> updateEdgeFindQuery;
    std::map<Edge,std::vector<std::pair<uint,uint>>> initEdgeFindQuery;
   
    std::vector<int>match;
    std::vector<std::vector<uint>> matchCandidate;
    std::vector<uint> indexCheckCNT;
    


    //combine part
    std::vector<std::vector<std::vector<uint>>> intersectionResult;
    std::vector<std::vector<int>> combineStack;
    std::vector<std::vector<uint>> headRecord;
    std::vector<uint> stackSize;
    std::vector<unFreezeStackType> type;
    std::vector<int> stackHead;
    bool print_init;

public:
    CSMPP(Graph& data_graph, Graph& query_graph, std::vector<Graph> & multiQueryGraph, uint max_num_results,
            bool print_prep, bool print_enum, bool homo, bool print_init);
    ~CSMPP() override{};

    void Preprocessing() override;
    void InitialMatching() override;

    void AddEdge(uint v1, uint v2, uint label) override;
    void RemoveEdge(uint v1, uint v2) override;
    void AddVertex(uint id, uint label) override;
    void RemoveVertex(uint id) override;
    
    void GetMemoryCost(size_t &num_edges, size_t &num_vertices) override;
    void TimePrint();
 
    

private:

    bool indexCheck(uint data_v, uint query_v, uint queryID);
    bool indexCheck(uint data_v, uint query_v, uint queryID, uint pos);
    void EdgeFindQueryInit();
    void initEdgeFindQueryInit();
    std::vector<std::pair<uint, uint>> UpdateEdgeFineQuery(uint v1, uint v2, uint edgeLabel, searchType type);
    void setMatchVertex(const std::vector<uint> & matchingIndex, const std::vector<int> & vertexs);
    void setMatchVertex(const std::vector<std::pair<uint,uint>> & matchingIndex, const std::vector<int> & vertexs);
    void unsetMatchVertex(const std::vector<uint> & matchingIndex);
    void unsetMatchVertex(const std::vector<std::pair<uint,uint>> & matchingIndex);
    void addMatchResult(uint queryIndex, uint edgeIndex, searchType type);
    #if defined(ADD_RESULT_DFS)
    void addMatchResult_DFS(uint queryIndex, uint edgeIndex, searchType type);
    
    void DFS_SEARCH_ALL(const std::vector<uint> & needToCombine, const std::vector<uint> & cacheMap, uint idx, searchType type);
    void DFS_SEARCH_ALL(const std::vector<std::vector<int>> & needToCombine, const std::vector<uint> & NoOverLeafWeight, const std::vector<uint> & isolatedIndex, uint idx, searchType type, size_t weight);
    #endif
    void printMatch();
    void searchInit(uint v1, uint v2, uint label, searchType type);
    void searchVertex(uint queryIndex, uint edgeIndex, searchType type, uint depth);
    void matchVertex(uint vertex, uint pos);
    void matchVertex(std::vector<uint> & Candidate, uint pos);
    
    void matchBatchVertices(uint index, std::vector<uint> & candidates);
    void unmatchBatchVertices(uint index);
    void popVertex(uint vertex, uint pos);
    void popVertex(uint pos);
    
    void printMatch(const std::vector<uint> & matchOrder);


    bool vertexPushCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex, uint elabel);
    bool getIntersection(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex, uint elabel);
    std::vector<uint> & getItersectionTop(int depth);
    void combineStackPopTail(int depth);
    void combineStackPopTail(int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight);
    void combinePushBack(uint vertex, uint vertexInWhichCandidateNeighbor, int depth);
    void combinePushBack(int vertex, uint vertexInWhichCandidateNeighbor, int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight);
    bool headChange(const std::vector<std::vector<uint>> & needToCombine, int depth);
    bool headChange(const std::vector<std::vector<int>> & needToCombine, int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight);

    void setVisitedPatch(const std::vector<int> & vertex);
    void setUnVisitedPatch(const std::vector<int> & vertex);

    static bool cmp(const std::pair<int,int> &p1,const std::pair<int,int> &p2){
        return p1.second < p2.second;
    }
};
