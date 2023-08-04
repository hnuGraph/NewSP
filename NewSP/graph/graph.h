#ifndef GRAPH_GRAPH
#define GRAPH_GRAPH

#include <queue>
#include <tuple>
#include <vector>
#include <map>
#include "utils/types.h"
#include "utils/utils.h"
#include "configuration/config.h"
#include "DecisionMakingSystem/DecisionMakingSystem.h"
#define LABLE_TYPE std::vector<std::pair<std::vector<uint>, std::vector<uint>>> 

class Edge{
private:
    uint v1;
    uint v2;
    uint v1Label;
    uint v2Label;
    uint eLabel;
    uint index;
    bool exist;
public:
    Edge(uint v1, uint v2, uint v1Label, uint v2Label, uint eLabel, uint index):v1(v1),v2(v2),v1Label(v1Label),v2Label(v2Label),eLabel(eLabel),index(index){
        this->exist = true;
    }
    Edge(uint v1Label, uint v2Label, uint eLabel):v1Label(v1Label),v2Label(v2Label),eLabel(eLabel){}
    bool GetExist(){return this->exist;}
    void EdgeDelete(){this->exist = false;}
    const uint GetV1() const{return this->v1;}
    const uint GetV2() const{return this->v2;}
    const uint GetV1Label() const{return this->v1Label;}
    const uint GetV2Label() const{return this->v2Label;}
    const uint GeteLabel() const{return this->eLabel;}
    uint GetIndex(){return this->index;}
    bool operator == (const Edge& edge)const {
        if((edge.v1Label == this->v1Label && edge.v2Label == this->v2Label && edge.eLabel == this->eLabel) ||
        (edge.v2Label == this->v1Label && edge.v1Label == this->v2Label && edge.eLabel == this->eLabel)){
            return true;
        }
        return false;
    }
    bool operator < (const Edge & edge) const{
        if(this->v1Label < edge.v1Label){
            return true;
        } else if (this->v1Label > edge.v1Label) {
            return false;
        }

        if(this->v2Label < edge.v2Label){
            return true;
        } else if (this->v2Label > edge.v2Label) {
            return false;
        }


        if(this->eLabel < edge.eLabel){
            return true;
        } else if (this->eLabel > edge.eLabel) {
            return false;
        }

        return false;
    }
};
class Graph
{
protected:
    uint edge_count_;
    uint vlabel_count_;
    uint elabel_count_;
    std::vector<std::vector<uint>> neighbors_;
    std::vector<std::vector<uint>> elabels_;
    std::vector<std::vector<uint>> matchOrder;
    std::vector<std::vector<uint>> matchCacheOrder;
    std::vector<std::vector<bool>> matchCacheSave;
    std::vector<std::vector<uint>> cacheMap;
    std::vector<std::vector<uint>> matchIsolatedOrder;
    std::vector<std::vector<vertexType>> matchVertexTypes;
    std::vector<std::vector<LRAndIndexCheckType>> DecisionList;
    std::vector<std::vector<std::vector<std::pair<uint,uint>>>> unfreezeRecord;
    std::vector<std::vector<std::vector<std::tuple<uint,uint,uint>>>> descList;//tuple <orderIndex, neighborLabel, elabel>
    std::vector<Edge> edge;
    std::vector<uint> mapping;
    std::vector<std::vector<uint>> freezeVertexNumAfter;
    std::vector<std::vector<uint>> isolatedVertexNumAfter;
    
public:
    std::queue<InsertUnit> updates_;
    std::vector<uint> vlabels_;// vertex's label
    bool matchOrderBuild;// build index flag
    LABLE_TYPE data_label; //data graph label check
#if GRAPH_TYPE == 0
    std::vector<int*>index;
#elif GRAPH_TYPE == 1
    std::vector<std::map<uint,int>> index;
#endif
    //make queryGraph use
    std::vector<uint>vertexLabel;
    std::vector<uint>vertexDataGraphID;
    std::map<uint, uint> isolatedVertexTimes;
    DMS decisionMakeSystem;

public:
    Graph()
    : edge_count_(0)
    , vlabel_count_(0)
    , elabel_count_(0)
    , neighbors_{}
    , elabels_{}
    , updates_{}
    , vlabels_{}
    , matchOrderBuild(false)
    {}
    Graph(bool matchOrderBuild):matchOrderBuild(matchOrderBuild){}

    virtual uint NumVertices() const { return vlabels_.size(); }
    virtual uint NumEdges() const { return edge_count_; }
    uint NumVLabels() const { return vlabel_count_; }
    uint NumELabels() const { return elabel_count_; }
    uint GetDiameter() const;

    void AddVertex(uint id, uint label);
    void RemoveVertex(uint id);
    void AddEdge(uint v1, uint v2, uint label);
    void RemoveEdge(uint v1, uint v2);

    uint GetVertexLabel(uint u) const;
    const std::vector<uint>& GetNeighbors(uint v) const;
    const std::vector<uint>& GetNeighborLabels(uint v) const;
    uint GetDegree(uint v) const;
    std::tuple<uint, uint, uint> GetEdgeLabel(uint v1, uint v2) const;

    void LoadFromFile(const std::string &path);
    void LoadFromFile(const std::string &prefixPath, std::vector<Graph> & queryGraph);
    void LoadUpdateStream(const std::string &path);
    void PrintMetaData() const;
    void PrintGraphNature() const;
    void indexAvg() const;

    const std::vector<uint> &  GetMatchOrder(uint index) const;
    void DescListInit(uint edgeIndex, std::vector<uint> & order, uint visitedVertex);
    std::vector<std::tuple<uint, uint, uint>> GetDescList(uint edgeIndex, uint depth);
    void MatchOrderInit();
    void EdgeMakeOrder(uint edgeIndex, uint v1, uint v2, std::vector<bool> & visit, std::vector<uint> & order);
    const std::vector<uint> & EdgeGetOrder(Edge edge) const;
    const std::vector<uint> & EdgeGetOrder(uint index) const;
    void matchOrderTypeSet();
    const vertexType getVertexType(uint edgeIndex, uint pos) const; 
    const uint getFreezeVertexNumAfter(uint edgeIndex, uint pos) const ;
    const uint getIsolatedVertexNumAfter(uint edgeIndex, uint pos) const;
    const std::vector<vertexType> & getVertexType(uint edgeIndex) const;
    void setVertexFree(uint edgeIndex, uint pos);
    const std::vector<std::pair<uint, uint>> & getUnfreezeList(uint edgeIndex, uint pos) const;
    void setVertexStatus(uint edgeIndex, const std::vector<std::pair<uint,uint>> & vertexs, vertexType type);
    const std::vector<uint>& getIsolatedVertexIndex(uint edgeIndex) const;
    const uint getDestListSize(uint egdeIndex, uint depth) const;
    const uint getUnfreezeListSize(uint egdeIndex, uint depth) const;

    uint getCacheMap(uint egdeIndex, uint VertexIndex);
    const std::vector<uint> & getCacheMap(uint egdeIndex) const;
    void initCacheOrder();
    void initDecision();
    LRAndIndexCheckType getDecision(uint EdgeIndex, uint VertexIndex);
    uint getCacheStatu(uint edgeIndex, uint depth);
    const std::vector<uint> & getCacheStatu(uint edgeIndex);
    bool getSaveStatu(uint edgeIndex, uint depth);
    void indexUpdate(uint v1, uint v2, uint label, bool op);//true insert; false delete
    void indexInit();
    void indexPrint();
    void MatchOrderAllPrint();
    void DescListPrint(uint edgeIndex);
    const uint getIndexValue(uint v, uint pos) const;

    Edge GetEdge(uint index);
    std::vector<Edge> GetEdge();

    bool mappingAdd(uint vertex, uint beginEdge);
    void mappingPopTail();
    uint GetMappingSize();

    void isolatedVertexTimesAdd(const std::vector<uint> & candidateVertexs);
    void isolatedVertexTimesMinus(const std::vector<uint> & candidateVertexs);


    static bool cmp(const std::pair<int,int> &p1,const std::pair<int,int> &p2){
        return p1.second > p2.second;
    }
};

#endif //GRAPH_GRAPH
