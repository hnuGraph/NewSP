#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <time.h>
#include <assert.h>
#include <bitset>

#include "utils/types.h"
#include "utils/globals.h"
#include "utils/utils.h"
#include "graph/graph.h"
#include "matching/CSMPP.h"
#include "configuration/config.h"


CSMPP::CSMPP(Graph& data_graph, 
        Graph& query_graph,
        std::vector<Graph> & multiQueryGraph, 
        uint max_num_results,
        bool print_prep, 
        bool print_enum, 
        bool homo,
        bool print_init
)
: matching(
    multiQueryGraph, query_graph, data_graph,  max_num_results, 
    print_prep, print_enum, homo)
{
    this->intersectionResult.resize(100);
    this->combineStack.resize(100);
    this->headRecord.resize(100);
    this->stackSize.resize(100);
    this->type.resize(100);
    this->stackHead.resize(100);
    this->print_init = print_init;
}
/**
 * @description: build matching order && make query-independnt index && cache index && adaptive dicision make
 * @return {*}
 */
void CSMPP::Preprocessing()
{
    this->data_.indexInit();
    for(auto & qGraph : this->queryVec)
    {
        qGraph.indexInit();
        qGraph.MatchOrderInit();
        qGraph.matchOrderTypeSet();
        qGraph.initDecision();
        qGraph.initCacheOrder();
    }
    this->EdgeFindQueryInit();
    if(this->print_init)
        this->initEdgeFindQueryInit();
}

/**
 * @description: query-independent cover index check
 * @param {uint} data_v
 * @param {uint} query_v
 * @param {uint} queryID
 * @return {*}
 */
bool CSMPP::indexCheck(uint data_v, uint query_v, uint queryID){
    const uint dataGraphVLabelMaxNum = this->data_.NumVLabels();
    const uint dataGraphELabelMaxNum = this->data_.NumELabels();
    const uint queryGraphVLabelMaxNum = this->queryVec[queryID].NumVLabels();
    const uint queryGraphELabelMaxNum = this->queryVec[queryID].NumELabels();

#if GRAPH_TYPE == 0
    const auto & dataGraph = this->data_.index[data_v];
    const auto & queryGraph = this->queryVec[queryID].index[query_v];
    for(int i = 0; i < queryGraphVLabelMaxNum; i++){
        if(i < dataGraphVLabelMaxNum){
            if(dataGraph[i] < queryGraph[i]){
                return false;
            }
        }
        else{
            if(queryGraph[i] > 0){
                return false;
            }
        }
    }
    for(int i = 0; i < queryGraphELabelMaxNum; i++){
        if(i < dataGraphELabelMaxNum){
            if(dataGraph[dataGraphVLabelMaxNum + i] < queryGraph[queryGraphVLabelMaxNum + i]){
                return false;
            }
        }
        else{
            if(queryGraph[queryGraphVLabelMaxNum + i] > 0){
                return false;
            }
        }
    }
    return true; 
    //for dense graph
#elif GRAPH_TYPE == 1
    auto & dataGraph = this->data_.index[data_v];
    auto & queryGraph = this->queryVec[queryID].index[query_v];
    auto iter = queryGraph.begin();
    while(iter != queryGraph.end())
    {
        auto dataIter = dataGraph.find(iter->first);
        if(dataIter == dataGraph.end() || dataIter->second < iter->second)
            return false;
        iter++;
    }
    return true;
#endif
}

/**
 * @description: query-independent cover index check
 * @param {uint} data_v
 * @param {uint} query_v
 * @param {uint} queryID
 * @param {uint} pos
 * @return {*}
 */
bool CSMPP::indexCheck(uint data_v, uint query_v, uint queryID, uint pos){
#if GRAPH_TYPE == 0
    this->indexCheckCNT[pos]++;
    const uint dataGraphVLabelMaxNum = this->data_.NumVLabels();
    const uint dataGraphELabelMaxNum = this->data_.NumELabels();
    const uint queryGraphVLabelMaxNum = this->queryVec[queryID].NumVLabels();
    const uint queryGraphELabelMaxNum = this->queryVec[queryID].NumELabels();
    const auto & dataGraph = this->data_.index[data_v];
    const auto & queryGraph = this->queryVec[queryID].index[query_v];
    for(int i = 0; i < queryGraphVLabelMaxNum; i++){
        if(i < dataGraphVLabelMaxNum){
            if(dataGraph[i] < queryGraph[i]){
                return false;
            }
        }
        else{
            if(queryGraph[i] > 0){
                return false;
            }
        }
    }
    for(int i = 0; i < queryGraphELabelMaxNum; i++){
        if(i < dataGraphELabelMaxNum){
            if(dataGraph[dataGraphVLabelMaxNum + i] < queryGraph[queryGraphVLabelMaxNum + i]){
                return false;
            }
        }
        else{
            if(queryGraph[queryGraphVLabelMaxNum + i] > 0){
                return false;
            }
        }
    }
    return true;
#endif
}

/**
 * @description: initial query subgraph matching process
 * @return {*}
 */
void CSMPP::InitialMatching()
{
    for(int i = 0; i < this->data_.NumEdges(); ++i){
        const auto & edge = this->data_.GetEdge(i);
        uint v1 = edge.GetV1();
        uint v2 = edge.GetV2();
        uint elabel = edge.GeteLabel();
        this->searchInit(v1, v2, elabel, init);
        if(reach_time_limit){
            break;
        }
    }
}

/**
 * @description: after a edge insert, update graph and begin search subgraph
 * @param {uint} v1
 * @param {uint} v2
 * @param {uint} label
 * @return {*}
 */
void CSMPP::AddEdge(uint v1, uint v2, uint label)
{
    this->data_.AddEdge(v1, v2, label);
    this->data_.indexUpdate(v1, v2, label, true);
    this->searchInit(v1, v2, label, pos);
}

/**
 * @description: after a edge remove, begin search subgraph and update graph
 * @param {uint} v1
 * @param {uint} v2
 * @return {*}
 */
void CSMPP::RemoveEdge(uint v1, uint v2)
{
    uint label = std::get<2>(this->data_.GetEdgeLabel(v1, v2));
    this->searchInit(v1, v2, label, neg);
    this->data_.RemoveEdge(v1, v2);
    this->data_.indexUpdate(v1, v2, label, false);
}

void CSMPP::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
    
    visited_.resize(id + 1, false);
}

void CSMPP::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void CSMPP::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;
}

/**
 * @description: build map for update edge to batch query graphs 
 * @return {*}
 */
void CSMPP::EdgeFindQueryInit(){
    for(uint i = 0; i < this->queryVec.size(); i++){
        Graph & query = this->queryVec[i];
        for(auto & edge : query.GetEdge()){
            uint label1 = edge.GetV1Label();
            uint label2 = edge.GetV2Label();
            if(label1 > label2)  {
                std::swap(label1, label2);
            }
            Edge triggerEdge(label1, label2, edge.GeteLabel());
            this->updateEdgeFindQuery[triggerEdge].emplace_back(std::make_pair(i, edge.GetIndex()));
        }
    }
}

/**
 * @description: build map for a vertex to batch query graphs
 * @return {*}
 */
void CSMPP::initEdgeFindQueryInit(){
    for(uint i = 0; i < this->queryVec.size(); i++){
        const Edge & edge = this->queryVec[i].GetEdge(0);
        uint label1 = edge.GetV1Label();
        uint label2 = edge.GetV2Label();
        if(label1 > label2){
            std::swap(label1, label2);
        }
        Edge triggerEdge(label1, label2, edge.GeteLabel());
        this->initEdgeFindQuery[triggerEdge].push_back(std::make_pair(i, 0));
    }
}

/**
 * @description: after graph update begin subgraph search
 * @param {uint} v1
 * @param {uint} v2
 * @param {uint} edgeLabel
 * @param {searchType} type
 * @return {*}
 */
std::vector<std::pair<uint, uint>> CSMPP::UpdateEdgeFineQuery(uint v1, uint v2, uint edgeLabel, searchType type){
    uint v1Label = this->data_.GetVertexLabel(v1);
    uint v2Label = this->data_.GetVertexLabel(v2);
    if(v1Label > v2Label){
        std::swap(v1Label, v2Label);
    }
    Edge edge(v1Label, v2Label, edgeLabel);
    if(type == pos || type == neg){
        if(this->updateEdgeFindQuery.find(edge) != this->updateEdgeFindQuery.end()){
            return this->updateEdgeFindQuery[edge];
        }
        else{
            return {};
        }
    }
    else{
        if(this->initEdgeFindQuery.find(edge) != this->initEdgeFindQuery.end()){
            return this->initEdgeFindQuery[edge];
        }
        else{
            return {};
        }
    }
}
/**
 * @description: begin subrgaph search
 * @param {uint} v1
 * @param {uint} v2
 * @param {uint} label
 * @param {searchType} type
 * @return {*}
 */
void CSMPP::searchInit(uint v1, uint v2, uint label, searchType type){
    //1. get a series queryTermLabel contain the update edge
    const auto &  querySeries = this->UpdateEdgeFineQuery(v1, v2, label, type);

    std::vector<std::pair<uint, uint>> queryCandidate;
    //2.index Check
    for(int index = 0; index < querySeries.size(); index++){
        uint queryGraph = querySeries[index].first;
        uint queryEdge = querySeries[index].second;
        Edge _queryEdge = this->queryVec[queryGraph].GetEdge(queryEdge);
        uint v1label = this->data_.GetVertexLabel(v1);
        uint v2label = this->data_.GetVertexLabel(v2);
        uint v1Query = _queryEdge.GetV1();
        uint v2Query = _queryEdge.GetV2();
        if(v1label == v2label){
            if((this->indexCheck(v1, v1Query, queryGraph) && this->indexCheck(v2, v2Query, queryGraph))||
            (this->indexCheck(v1, v2Query, queryGraph) && this->indexCheck(v2, v1Query, queryGraph))){
                queryCandidate.emplace_back(querySeries[index]);
            }
        }
        else{
            if(v1label == _queryEdge.GetV1Label()){
                if(this->indexCheck(v1, v1Query, queryGraph) &&
                    this->indexCheck(v2, v2Query, queryGraph)
                ) {
                    queryCandidate.emplace_back(querySeries[index]);
                }
            }
            else{
                if(this->indexCheck(v2, v1Query, queryGraph) &&
                    this->indexCheck(v1, v2Query, queryGraph)
                ) {
                    queryCandidate.emplace_back(querySeries[index]);
                }
            }
        }
    }
    uint64_t temp = 0;
    //3. begin seach
    for(auto & item : queryCandidate){
        auto  _edge = this->queryVec[item.first].GetEdge(item.second);
        uint _edgeV1label = _edge.GetV1Label();
        uint _edgeV2label = _edge.GetV2Label();
        const auto & matchOrder = this->queryVec[item.first].GetMatchOrder(item.second);
        if(_edgeV1label != _edgeV2label){
            if(this->data_.GetVertexLabel(v1) != this->queryVec[item.first].GetVertexLabel(matchOrder[0])){
                std::swap(v1, v2);
            }
            this->match.resize(this->queryVec[item.first].NumVertices(), UNMATCHED);
            this->matchCandidate.resize(this->queryVec[item.first].NumVertices());
            this->matchVertex(v1, 0);
            this->matchVertex(v2, 1); 
            this->queryVec[item.first].isolatedVertexTimes.clear();
            this->searchVertex(item.first, item.second, type, 2);

            this->popVertex(v2, 1);
            this->popVertex(v1, 0);
        }
        else{
            for(int i = 0; i < 2; i++){
                //index check
                if(this->indexCheck(v1, matchOrder[0], item.first)){
                    this->match.resize(this->queryVec[item.first].NumVertices(), UNMATCHED);
                    this->matchCandidate.resize(this->queryVec[item.first].NumVertices());
                    this->matchVertex(v1, 0);
                    this->matchVertex(v2, 1);
                    this->queryVec[item.first].isolatedVertexTimes.clear();
                    this->searchVertex(item.first, item.second, type, 2);
                    this->popVertex(v2, 1);
                    this->popVertex(v1, 0);
                }
                std::swap(v1,v2);// round 2 need
            }
        }
        if(reach_time_limit) return;
    }
}

/**
 * @description: search subgraph match step by step
 * @param {uint} queryIndex
 * @param {uint} edgeIndex
 * @param {searchType} type
 * @param {uint} depth
 * @return {*}
 */
void CSMPP::searchVertex(uint queryIndex, uint edgeIndex, searchType type, uint depth){
    auto & queryGraph = this->queryVec[queryIndex];
    const uint cacheStatu = queryGraph.getCacheStatu(edgeIndex, depth);
    const auto & matchOrder = queryGraph.GetMatchOrder(edgeIndex);
    const uint queryVexter = matchOrder[depth];
    vertexType currentSearchVertexType = queryGraph.getVertexType(edgeIndex, depth);
    const auto & desItem = queryGraph.GetDescList(edgeIndex, depth);//<v1Index,v1Label,eLabel>
    const auto & freezeIndex = queryGraph.getUnfreezeList(edgeIndex, depth);
    useCacheStatus finish = notUse;
    //use cache?
    if(depth != cacheStatu)
    {

        if(desItem.size() == 0 && freezeIndex.size() == 0)
        {
            // std::cout << "cache use finish" << std::endl;
            finish = useAndFinsh;
            if(currentSearchVertexType == freeVertex)
            {
                for(auto & i : this->matchCandidate[cacheStatu])
                {
                    if(this->visited_[i])continue;
                    this->matchVertex(i,depth);
                    this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
                    this->popVertex(i, depth);
                }
            }
            else{
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesAdd(this->matchCandidate[cacheStatu]);
                }
                if(depth == queryGraph.NumVertices() - 1)
                {
                #if defined(ADD_RESULT_DFS)
                    this->addMatchResult_DFS(queryIndex, edgeIndex, type);
                #else
                    this->addMatchResult(queryIndex, edgeIndex, type);
                #endif
                }
                else
                {
                    this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
                }
                if(currentSearchVertexType == isolatedVertex){
                    this->queryVec[queryIndex].isolatedVertexTimesMinus(this->matchCandidate[cacheStatu]);
                }
            }
        }
        else
        {
            finish = useNotFinsh;
        }
    }
    
    uint vertexLabel = queryGraph.GetVertexLabel(queryVexter);

    //1.intersection worst case
    //1.1 find min
    if(finish != useAndFinsh)
    {
    std::vector<uint> candidate;
    LRAndIndexCheckType decision = queryGraph.getDecision(edgeIndex, depth); // adaptive label distribute check
    if(desItem.size() != 0){
        uint min_u_index = INT_MAX;
        uint min_u_neighborSize = INT_MAX;
        for(int i = 0; i < desItem.size(); i++){
            const auto & item = desItem[i];
            uint v = this->match[std::get<0>(item)];
            uint Size = this->data_.getIndexValue(v, vertexLabel);
            if(Size < min_u_neighborSize){
                min_u_index = i;
                min_u_neighborSize = Size;
            }
            else if(Size == min_u_neighborSize){
                uint preID = this->match[std::get<0>(desItem[min_u_index])];
                if(this->data_.GetNeighbors(preID).size() > this->data_.GetNeighbors(v).size()){
                    min_u_index = i;
                    min_u_neighborSize = Size;
                }
            }
        }
        //1.use cache? 2. cache is the min set? 
        if(finish == useNotFinsh)
        {
            if(min_u_neighborSize < this->matchCandidate[cacheStatu].size())
            {
                uint min_u = this->match[std::get<0>(desItem[min_u_index])];
                uint min_u_elabel = std::get<2>(desItem[min_u_index]);//elabel
                const auto& MinNeighbor = this->data_.GetNeighbors(min_u);
                const auto& q_nbr_labels = this->data_.GetNeighborLabels(min_u);
                for(int i = 0; i < MinNeighbor.size(); i++){
                    const uint v = MinNeighbor[i];
                    //1. check label
                    if(
                        this->data_.GetVertexLabel(v) != vertexLabel ||
                        q_nbr_labels[i] != min_u_elabel
                    )continue;
                    //2. check visit
                    if(this->visited_[v] == true && !homomorphism_ ) continue;
                    //3.check if joinable
                    bool joinable = true;
                    for(int k = 0; k < desItem.size(); k++){
                        if(k == min_u_index) continue;
                        const uint data_V = this->match[std::get<0>(desItem[k])];
                        const uint elabel = std::get<2>(desItem[k]);
                        const auto & dataVNeighbor = this->data_.GetNeighbors(data_V);
                        auto it = std::lower_bound(dataVNeighbor.begin(), dataVNeighbor.end(), v);
                        if(it == dataVNeighbor.end() || 
                        *it != v
                        ){
                            joinable = false;
                            break;
                        }
                        else
                        {
                            uint dis = std::distance(dataVNeighbor.begin(), it);
                            if(this->data_.GetNeighborLabels(data_V)[dis] != elabel)
                            {
                                joinable = false;
                                break;
                            }
                        }
                    }
                    auto it = std::lower_bound(this->matchCandidate[cacheStatu].begin(), this->matchCandidate[cacheStatu].end(), v);
                    if(it == this->matchCandidate[cacheStatu].end() || 
                    *it != v
                    ){
                        joinable = false;
                    }
                    if(!joinable)continue;
                    if(decision == Part1Check && !this->indexCheck(v, queryVexter, queryIndex))continue;
                    candidate.emplace_back(v);
                }                
            }
            else
            {
                for(int i = 0; i < this->matchCandidate[cacheStatu].size(); i++){
                    const uint v = this->matchCandidate[cacheStatu][i];
                    //1. check visit
                    if(this->visited_[v] == true && !homomorphism_ ) continue;
                    //2.check if joinable
                    bool joinable = true;
                    for(int k = 0; k < desItem.size(); k++){
                        const uint data_V = this->match[std::get<0>(desItem[k])];
                        const uint elabel = std::get<2>(desItem[k]);
                        const auto & dataVNeighbor = this->data_.GetNeighbors(data_V);
                        auto it = std::lower_bound(dataVNeighbor.begin(), dataVNeighbor.end(), v);
                        if(it == dataVNeighbor.end() || 
                        *it != v
                        ){
                            joinable = false;
                            break;
                        }
                        else
                        {
                            uint dis = std::distance(dataVNeighbor.begin(), it);
                            if(this->data_.GetNeighborLabels(data_V)[dis] != elabel)
                            {
                                joinable = false;
                                break;
                            }
                        }
                    }
                    if(!joinable)continue;
                    if(decision == Part1Check && !this->indexCheck(v, queryVexter, queryIndex))continue;
                    candidate.emplace_back(v);
                }                
            }
            finish = useAndFinsh;
        }
        else
        {
            uint min_u = this->match[std::get<0>(desItem[min_u_index])];
            uint min_u_elabel = std::get<2>(desItem[min_u_index]);//elabel
            const auto& MinNeighbor = this->data_.GetNeighbors(min_u);
            const auto& q_nbr_labels = this->data_.GetNeighborLabels(min_u);
            for(int i = 0; i < MinNeighbor.size(); i++){
                const uint v = MinNeighbor[i];
                //1. check label
                if(
                    this->data_.GetVertexLabel(v) != vertexLabel ||
                    q_nbr_labels[i] != min_u_elabel
                )continue;
                //2. check visit
                if(this->visited_[v] == true && !homomorphism_ ) continue;
                //3.check if joinable
                bool joinable = true;
                for(int k = 0; k < desItem.size(); k++){
                    if(k == min_u_index) continue;
                    const uint data_V = this->match[std::get<0>(desItem[k])];
                    const uint elabel = std::get<2>(desItem[k]);
                    const auto & dataVNeighbor = this->data_.GetNeighbors(data_V);
                    auto it = std::lower_bound(dataVNeighbor.begin(), dataVNeighbor.end(), v);
                    if(it == dataVNeighbor.end() || 
                    *it != v
                    ){
                        joinable = false;
                        break;
                    }
                    else
                    {
                        uint dis = std::distance(dataVNeighbor.begin(), it);
                        if(this->data_.GetNeighborLabels(data_V)[dis] != elabel)
                        {
                            joinable = false;
                            break;
                        }
                    }
                }
                if(!joinable)continue;
                if(decision == Part1Check && !this->indexCheck(v, queryVexter, queryIndex))continue;
                candidate.emplace_back(v);
            }
        }
        if(candidate.size() == 0)
            return;
    }
    //expansion virtual node
    if(!freezeIndex.empty()){
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();  
        std::vector<uint> elabels;
        std::vector<std::vector<uint>> needToCombine;
        needToCombine.resize(freezeIndex.size());
        for(int i = 0; i < freezeIndex.size(); i++){
            needToCombine[i].reserve(this->matchCandidate[freezeIndex[i].first].size());
            if(this->queryVec[queryIndex].getDestListSize(edgeIndex, freezeIndex[i].first) == 0){
                for(auto item : this->matchCandidate[queryGraph.getCacheStatu(edgeIndex, freezeIndex[i].first)])
                {
                    if(this->indexCheck(item, matchOrder[freezeIndex[i].first], queryIndex))
                    {
                        needToCombine[i].emplace_back(item);
                    }
                }
                if(needToCombine[i].size() == 0){
                    return;//all kill
                }
            }
            else
            {
                needToCombine[i] = this->matchCandidate[freezeIndex[i].first];
            }
            this->headRecord[depth].emplace_back(needToCombine[i].size());
            elabels.emplace_back(freezeIndex[i].second);
        }
        //combine part
        //init
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = runningStack;
        //combine
        //int stackResult = 0;
        while(this->headRecord[depth][0] >= 0){
            //1.push back stack
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth);
                    if(pending == false){
                        return;
                    }
                }
                int replaceIndex = this->stackHead[depth];// need to add
                uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(finish != useNotFinsh)
                {
                    if(vertexPushCheck(currentVertex, vertexLabel, candidate, depth, queryIndex, queryVexter, elabels[replaceIndex]) == false){
                        this->headRecord[depth][replaceIndex]--;
                    }
                    else{
                        this->combinePushBack(currentVertex, replaceIndex, depth);
                    }
                }
                else
                {
                    if(vertexPushCheck(currentVertex, vertexLabel, this->matchCandidate[cacheStatu], depth, queryIndex, queryVexter, elabels[replaceIndex]) == false){
                        this->headRecord[depth][replaceIndex]--;
                    }
                    else{
                        this->combinePushBack(currentVertex, replaceIndex, depth);
                    }
                }
            }
            //2.search
            this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freeVertex);
            this->setMatchVertex(freezeIndex, this->combineStack[depth]);
            if(decision == Part2Check)
            {
                std::vector<uint> Part2Candidate;
                for(auto & dataV : this->getItersectionTop(depth)){
                    if(this->indexCheck(dataV, queryVexter, queryIndex)){
                        Part2Candidate.emplace_back(dataV);
                    }
                }
                if(Part2Candidate.size() != 0){
                    if(currentSearchVertexType == freeVertex){
                        for(auto & dataV : Part2Candidate){
                            this->matchVertex(dataV, depth);
                            this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
                            this->popVertex(dataV, depth);
                        }
                    }
                    else{    
                        this->matchVertex(Part2Candidate, depth);
                        if(currentSearchVertexType == isolatedVertex){
                            this->queryVec[queryIndex].isolatedVertexTimesAdd(this->matchCandidate[depth]);
                        }
                        if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                        #if defined(ADD_RESULT_DFS)
                            this->addMatchResult_DFS(queryIndex, edgeIndex, type);
                        #else
                            this->addMatchResult(queryIndex, edgeIndex, type);
                        #endif
                        }
                        else{
                            this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
                        }
                        if(currentSearchVertexType == isolatedVertex){
                            this->queryVec[queryIndex].isolatedVertexTimesMinus(this->matchCandidate[depth]);
                        }
                        this->popVertex(depth);
                    }
                }
            }
            else{
                if(currentSearchVertexType == freeVertex){
                    for(auto & dataV : this->getItersectionTop(depth)){
                        this->matchVertex(dataV, depth);
                        this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
                        this->popVertex(dataV, depth);
                    }
                }
                else{
                    this->matchVertex(this->getItersectionTop(depth), depth);
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesAdd(this->matchCandidate[depth]);
                    }
                    if(depth == this->queryVec[queryIndex].NumVertices() - 1){
                    #if defined(ADD_RESULT_DFS)
                        this->addMatchResult_DFS(queryIndex, edgeIndex, type);
                    #else
                        this->addMatchResult(queryIndex, edgeIndex, type);
                    #endif
                    }
                    else{
                        this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
                    }
                    if(currentSearchVertexType == isolatedVertex){
                        this->queryVec[queryIndex].isolatedVertexTimesMinus(this->matchCandidate[depth]);
                    }
                    this->popVertex(depth);
                }
            }
            this->unsetMatchVertex(freezeIndex);
            this->queryVec[queryIndex].setVertexStatus(edgeIndex, freezeIndex, freezeVertex);
            this->combineStackPopTail(depth);
        }
    }
    else{
        if(currentSearchVertexType == freeVertex){
            for(auto & dataV : candidate){
                this->matchVertex(dataV, depth);
                this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
                this->popVertex(dataV, depth);
            }
        }
        else{
            this->matchVertex(candidate, depth);
            if(currentSearchVertexType == isolatedVertex){
                this->queryVec[queryIndex].isolatedVertexTimesAdd(this->matchCandidate[depth]);
            }
            if(depth == this->queryVec[queryIndex].NumVertices() - 1){
            #if defined(ADD_RESULT_DFS)
                this->addMatchResult_DFS(queryIndex, edgeIndex, type);
            #else
                this->addMatchResult(queryIndex, edgeIndex, type);
            #endif
            }
            else{
                this->searchVertex(queryIndex, edgeIndex, type, depth + 1);
            }
            if(currentSearchVertexType == isolatedVertex){
                this->queryVec[queryIndex].isolatedVertexTimesMinus(this->matchCandidate[depth]);
            }
            this->popVertex(depth);
        }
    }
    }
}

/**
 * @description: make batch query vertices match
 * @param {vector<uint>} &
 * @return {*}
 */
void CSMPP::setMatchVertex(const std::vector<uint> & matchingIndex, const std::vector<int>& vertexs){
    for(int i = 0; i < matchingIndex.size(); i++){
        //assert(this->visited_[vertexs[i]] == false);
        this->match[matchingIndex[i]] = vertexs[i];
        this->visited_[vertexs[i]] = true;
    }
}

/**
 * @description: make batch query vertices match
 * @param {vector<uint>} &
 * @return {*}
 */
void CSMPP::setMatchVertex(const std::vector<std::pair<uint,uint>> & matchingIndex, const std::vector<int> & vertexs)
{
    for(int i = 0; i < matchingIndex.size(); i++){
        //assert(this->visited_[vertexs[i]] == false);
        this->match[matchingIndex[i].first] = vertexs[i];
    }
}

/**
 * @description: make batch query vertices unmatch
 * @param {vector<uint>} &
 * @return {*}
 */
void CSMPP::unsetMatchVertex(const std::vector<uint> & matchingIndex){
    for(int i = 0; i < matchingIndex.size(); i++){
        this->visited_[this->match[matchingIndex[i]]] = false;
        this->match[matchingIndex[i]] = -1;
    }
}

/**
 * @description: make batch query vertices unmatch
 * @param {vector<uint>} &
 * @return {*}
 */
void CSMPP::unsetMatchVertex(const std::vector<std::pair<uint,uint>> & matchingIndex)
{
    for(int i = 0; i < matchingIndex.size(); i++){
        this->match[matchingIndex[i].first] = -1;
    }
}

/**
 * @description: final make multi expansion
 * @param {uint} queryIndex
 * @param {uint} edgeIndex
 * @param {searchType} type
 * @return {*}
 */
void CSMPP::addMatchResult(uint queryIndex, uint edgeIndex, searchType type){
    const auto & isolateVertexIndex = this->queryVec[queryIndex].getIsolatedVertexIndex(edgeIndex);
    const auto & cacheMap = this->queryVec[queryIndex].getCacheMap(edgeIndex);
    if(print_enumeration_results_){
        std::vector<std::vector<uint>> needToCombine;
        for(int i = 0; i < isolateVertexIndex.size(); i++){
            needToCombine.emplace_back(this->matchCandidate[cacheMap[isolateVertexIndex[i]]]);
        }
        //combine part
        //init
        uint depth = this->queryVec[queryIndex].NumVertices();
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();
        for(int i = 0; i < needToCombine.size(); i++){
            headRecord[depth].push_back(needToCombine[i].size());
        }
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = finalStack;
        //combine
        while(this->headRecord[depth][0] >= 0){
            //1.push back stack
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth);
                    if(pending == false){
                        return;
                    }
                }
                int replaceIndex = this->stackHead[depth];// need to add
                uint currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(this->visited_[currentVertex] == true){
                    this->headRecord[depth][replaceIndex]--;
                }
                else{
                    this->combinePushBack(currentVertex, replaceIndex, depth);
                }
            }
            if(type == pos){
                num_positive_results_ ++;
            }
            else{
                if(type == neg){
                    num_negative_results_++;
                }
                else{
                    this->num_initial_results_++;
                }
            }
            this->setUnVisitedPatch(this->combineStack[depth]);
            this->setMatchVertex(isolateVertexIndex, this->combineStack[depth]);
            this->printMatch(this->queryVec[queryIndex].GetMatchOrder(edgeIndex));
            this->unsetMatchVertex(isolateVertexIndex);
            this->setVisitedPatch(this->combineStack[depth]);
            this->combineStackPopTail(depth);
        }
    }
    else{
        uint depth = this->queryVec[queryIndex].NumVertices();
        auto & isolatedVertexMap = this->queryVec[queryIndex].isolatedVertexTimes;
        std::vector<std::vector<int>> needToCombineV1;
        std::vector<uint> NoOverLeafWeight;
        bool allSame = true;
        for(int i = 0; i < isolateVertexIndex.size(); i++){
            const auto & I_isolateVertexCandidate = this->matchCandidate[cacheMap[isolateVertexIndex[i]]];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(this->visited_[vertex] == true){
                    continue;
                }
                I_needToCombine.push_back(vertex);
            }
            if(I_needToCombine.empty()){
                return;
            }
            needToCombineV1.push_back(I_needToCombine);
            if(allSame && i >= 1){
                if(needToCombineV1[i - 1].size() != needToCombineV1[i].size()){
                    allSame = false;
                }
            }
        }
        if(needToCombineV1.size() == 1){
            if(type == pos){
                num_positive_results_ += needToCombineV1[0].size();
            }
            else{
                if(type == neg){
                    num_negative_results_ += needToCombineV1[0].size();
                }
                else{
                    num_initial_results_ += needToCombineV1[0].size();
                }
            }
            
            return;
        }
        if(allSame == true){
            const auto & firstItem = needToCombineV1[0];
            for(int k = 1; k < needToCombineV1.size(); k++){
                const auto & kSelf = needToCombineV1[k];
                if(firstItem != kSelf){
                    allSame = false;
                }
            }
            if(allSame == true){
                size_t Matchresult = 1;
                for(int i = 0; i < needToCombineV1.size(); i++){
                    if(firstItem.size() - i <= 0)
                    {
                        return;
                    }
                    Matchresult *= (firstItem.size() - i);
                }
                if(type == pos){
                    num_positive_results_ += Matchresult;
                }
                else{
                    if(type == neg){
                        num_negative_results_ += Matchresult;
                    }
                    else{
                        num_initial_results_ += Matchresult;
                    }
                }
                
                return;
            }
        }
        std::vector<std::vector<int>> needToCombine;
        for(int i = 0; i < needToCombineV1.size(); i++){
            const auto & I_isolateVertexCandidate = needToCombineV1[i];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            int I_NoOverLeafWeight = 0;
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(isolatedVertexMap[vertex] > 1){
                    I_needToCombine.push_back(vertex);
                }
                else{
                    I_NoOverLeafWeight++;
                }
            }
            NoOverLeafWeight.push_back(I_NoOverLeafWeight);
            if(I_NoOverLeafWeight > 0){
                I_needToCombine.push_back(-1);
            }
            needToCombine.push_back(I_needToCombine);
        }

        //combine part
        //init
        this->headRecord[depth].clear();
        this->intersectionResult[depth].clear();
        this->combineStack[depth].clear();
        for(int i = 0; i < needToCombine.size(); i++){
            headRecord[depth].push_back(needToCombine[i].size());
        }
        this->stackSize[depth] = needToCombine.size();
        this->stackHead[depth] = 0;
        this->combineStack[depth].resize(this->stackSize[depth]);
        this->intersectionResult[depth].resize(this->stackSize[depth]);
        this->type[depth] = finalStack;
        //combine
        size_t totalMatch = 1;
        int currentVertex;
        int replaceIndex;
        while(this->headRecord[depth][0] >= 0){
            //1.push back stack
            while(this->stackHead[depth] < this->stackSize[depth]){
                if(this->headRecord[depth][this->stackHead[depth]] == 0){
                    bool pending = this->headChange(needToCombine, depth, totalMatch, NoOverLeafWeight);
                    if(pending == false){
                        return;
                    }
                }
                replaceIndex = this->stackHead[depth];// need to add
                currentVertex = needToCombine[replaceIndex][this->headRecord[depth][replaceIndex] - 1];
                if(currentVertex != -1 && this->visited_[currentVertex] == true){
                    this->headRecord[depth][replaceIndex]--;
                }
                else{
                    this->combinePushBack(currentVertex, replaceIndex, depth, totalMatch, NoOverLeafWeight);
                }
            }
            if(type == pos){
                num_positive_results_ += totalMatch;
            }
            else{
                if(type == neg){
                    num_negative_results_ += totalMatch;
                }
                else{
                    num_initial_results_ += totalMatch;
                }
                
            }
            
            this->combineStackPopTail(depth, totalMatch, NoOverLeafWeight);
        }
    }
}

/**
 * @description: show match result
 * @param {vector<uint>} &
 * @return {*}
 */
void CSMPP::printMatch(const std::vector<uint> & matchOrder){
    std::vector<int>matchCopy(match);
    std::cout << "vertex list" << std::endl;
    for(int i = 0; i < this->match.size(); i++){
        std::cout << this->match[i] << " ";
    }
    std::cout << std::endl;
}

/**
 * @description: match queryVertex
 * @param {uint} vertex
 * @param {uint} pos
 * @return {*}
 */
void CSMPP::matchVertex(uint vertex, uint pos){
    this->match[pos] = vertex;
    this->visited_[vertex] = true;
}

/**
 * @description: match batch of query vertices
 * @param {uint} vertex
 * @param {uint} pos
 * @return {*}
 */
void CSMPP::matchVertex(std::vector<uint> & Candidate, uint pos){
    this->match[pos] = -1;
    this->matchCandidate[pos] = std::move(Candidate);
}

/**
 * @description: unmatch query vertex
 * @param {uint} vertex
 * @param {uint} pos
 * @return {*}
 */
void CSMPP::popVertex(uint vertex, uint pos){
    this->match[pos] = UNMATCHED;
    this->visited_[vertex] = false;
}

/**
 * @description: unmatch query vertices
 * @param {uint} pos
 * @return {*}
 */
void CSMPP::popVertex(uint pos){
    this->match[pos] = UNMATCHED;
    this->matchCandidate[pos].clear();
}

/**
 * @description: check a combine of multi query whether can expansion
 * @param {uint} vertex
 * @param {uint} queryVertexLabel
 * @param {vector<uint>} &
 * @param {int} depth
 * @param {uint} queryIndex
 * @param {uint} queryVertex
 * @param {uint} elabel
 * @return {*}
 */
bool CSMPP::vertexPushCheck(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex, uint elabel){
    if(this->visited_[vertex] == true){
        return false;
    }
    if(this->type[depth] == runningStack){
        if(this->getIntersection(vertex, queryVertexLabel, needToIntersection, depth, queryIndex, queryVertex, elabel) == false){
            return false;
        }
    }
    return true;
}

/**
 * @description: get neighbor intersection result
 * @param {uint} vertex
 * @param {uint} queryVertexLabel
 * @param {vector<uint>} &
 * @param {int} depth
 * @param {uint} queryIndex
 * @param {uint} queryVertex
 * @param {uint} elabel
 * @return {*}
 */
bool CSMPP::getIntersection(uint vertex, uint queryVertexLabel, const std::vector<uint> & needToIntersection, int depth, uint queryIndex, uint queryVertex, uint elabel){
    std::vector<uint>result;
    const auto & vNeighbors = this->data_.GetNeighbors(vertex);
    const auto & vNeighborsLabels = this->data_.GetNeighborLabels(vertex);
    if(!needToIntersection.empty()){
        for(auto & item : needToIntersection)
        {
            auto it = std::lower_bound(vNeighbors.begin(), vNeighbors.end(), item);
            if(it != vNeighbors.end() && *it == item)
            {
                uint dis = std::distance(vNeighbors.begin(), it);
                if(vNeighborsLabels[dis] == elabel)
                {
                    result.emplace_back(item);
                }
            }
        }
    }
    else{
        for(int i = 0; i < vNeighbors.size(); i++){
            if(this->visited_[vNeighbors[i]] == false 
            &&this->data_.GetVertexLabel(vNeighbors[i]) == queryVertexLabel 
            &&vNeighborsLabels[i] == elabel
            ){
                result.emplace_back(vNeighbors[i]);
            }
        }
    }
    if(result.empty()){//intersecion result is null
        return false;
    }
    if(this->stackHead[depth] == 0){
        this->intersectionResult[depth][this->stackHead[depth]] = result;
        return true;
    }
    std::vector<uint> finalResult;
    const std::vector<uint> & preVec = this->getItersectionTop(depth);
    finalResult.reserve(preVec.size());
    std::set_intersection(preVec.begin(), preVec.end(), result.begin(), result.end(), std::back_inserter(finalResult));
    if(finalResult.empty()){
        return false;
    }
    this->intersectionResult[depth][this->stackHead[depth]] = std::move(finalResult);
    return true;
}

/**
 * @description: get intersecion result
 * @param {int} depth
 * @return {*}
 */
std::vector<uint> & CSMPP::getItersectionTop(int depth){
    return this->intersectionResult[depth][this->stackHead[depth] - 1];
}

/**
 * @description: pop out a match vertex
 * @param {int} depth
 * @return {*}
 */
void CSMPP::combineStackPopTail(int depth){
    uint vertex = this->combineStack[depth][this->stackHead[depth] - 1];
    this->stackHead[depth]--;
    this->visited_[vertex] = false;
}

/**
 * @description: pop out a match vertex
 * @param {int} depth
 * @return {*}
 */
void CSMPP::combineStackPopTail(int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight){
    uint vertex = this->combineStack[depth][this->stackHead[depth] - 1];
    
    if(vertex != -1){
        this->visited_[vertex] = false;
    }
    else{
        totalMatch /= NoOverLeftWeight[this->stackHead[depth] - 1];
    }
    this->stackHead[depth]--;
}

/**
 * @description: push back a match vertex
 * @param {uint} vertex
 * @param {uint} vertexInWhichCandidateNeighbor
 * @param {int} depth
 * @return {*}
 */
void CSMPP::combinePushBack(uint vertex, uint vertexInWhichCandidateNeighbor, int depth){
    this->visited_[vertex] = true;
    this->combineStack[depth][this->stackHead[depth]] = vertex;
    this->stackHead[depth]++;
    this->headRecord[depth][vertexInWhichCandidateNeighbor]--;
}

/**
 * @description: push back a match vertex
 * @param {uint} vertex
 * @param {uint} vertexInWhichCandidateNeighbor
 * @param {int} depth
 * @return {*}
 */
void CSMPP::combinePushBack(int vertex, uint vertexInWhichCandidateNeighbor, int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight){
    if(vertex != -1){
        this->visited_[vertex] = true;
    }
    else{
        totalMatch *= NoOverLeftWeight[vertexInWhichCandidateNeighbor];
    }
    this->combineStack[depth][this->stackHead[depth]] = vertex;
    this->stackHead[depth]++;
    this->headRecord[depth][vertexInWhichCandidateNeighbor]--;
}

/**
 * @description: reset combine query vertex
 * @param {vector<std::vector<uint>>} &
 * @param {int} depth
 * @return {*}
 */
bool CSMPP::headChange(const std::vector<std::vector<uint>> & needToCombine, int depth){
    int iHead = this->stackHead[depth];
    while(this->headRecord[depth][iHead] <= 0){
        if(iHead <= 0){
            return false;
        }
        this->headRecord[depth][iHead] = needToCombine[iHead].size();
        iHead--;
        this->combineStackPopTail(depth);
    }
    return true;
}

/**
 * @description: reset combine query vertex
 * @param {vector<std::vector<uint>>} &
 * @param {int} depth
 * @return {*}
 */
bool CSMPP::headChange(const std::vector<std::vector<int>> & needToCombine, int depth, size_t & totalMatch, const std::vector<uint>& NoOverLeftWeight){
    int iHead = this->stackHead[depth];
    while(this->headRecord[depth][iHead] <= 0){
        if(iHead <= 0){
            return false;
        }
        this->headRecord[depth][iHead] = needToCombine[iHead].size();
        iHead--;
        this->combineStackPopTail(depth, totalMatch, NoOverLeftWeight);
    }
    return true;
}

/**
 * @description: visit batch data vertices
 * @param {vector<int>} &
 * @return {*}
 */
void CSMPP::setVisitedPatch(const std::vector<int> & vertex){
    for(int i = 0; i < vertex.size(); i++){
        this->visited_[vertex[i]] = true;
    }
}

/**
 * @description: unvisited batch data vertices
 * @param {vector<int>} &
 * @return {*}
 */
void CSMPP::setUnVisitedPatch(const std::vector<int> & vertex){
    for(int i = 0; i < vertex.size(); i++){
        this->visited_[vertex[i]] = false;
    }
}
#if defined(ADD_RESULT_DFS)
/**
 * @description: add result using DFS
 * @param {uint} queryIndex
 * @param {uint} edgeIndex
 * @param {searchType} type
 * @return {*}
 */
void CSMPP::addMatchResult_DFS(uint queryIndex, uint edgeIndex, searchType type)
{
    const auto & isolateVertexIndex = this->queryVec[queryIndex].getIsolatedVertexIndex(edgeIndex);
    const auto & cacheMap = this->queryVec[queryIndex].getCacheMap(edgeIndex);

    if(print_enumeration_results_)
    {
        DFS_SEARCH_ALL(isolateVertexIndex, cacheMap, 0, type);
    }
    else
    {
        uint depth = this->queryVec[queryIndex].NumVertices();
        auto & isolatedVertexMap = this->queryVec[queryIndex].isolatedVertexTimes;
        std::vector<std::vector<int>> needToCombineV1;
        std::vector<uint> NoOverLeafWeight;
        bool allSame = true;
        for(int i = 0; i < isolateVertexIndex.size(); i++){
            const auto & I_isolateVertexCandidate = this->matchCandidate[cacheMap[isolateVertexIndex[i]]];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(this->visited_[vertex] == true){
                    continue;
                }
                I_needToCombine.push_back(vertex);
            }
            if(I_needToCombine.empty()){
                return;
            }
            needToCombineV1.push_back(I_needToCombine);
            if(allSame && i >= 1){
                if(needToCombineV1[i - 1].size() != needToCombineV1[i].size()){
                    allSame = false;
                }
            }
        }
        if(needToCombineV1.size() == 1){
            if(type == pos){
                num_positive_results_ += needToCombineV1[0].size();
            }
            else{
                if(type == neg){
                    num_negative_results_ += needToCombineV1[0].size();
                }
                else{
                    num_initial_results_ += needToCombineV1[0].size();
                }
            }
            
            return;
        }
        if(allSame == true){
            const auto & firstItem = needToCombineV1[0];
            for(int k = 1; k < needToCombineV1.size(); k++){
                const auto & kSelf = needToCombineV1[k];
                if(firstItem != kSelf){
                    allSame = false;
                    break;
                }
            }
            if(allSame == true){
                //std::cout << "allSame" << std::endl;
                size_t Matchresult = 1;
                for(int i = 0; i < needToCombineV1.size(); i++){
                    if(firstItem.size() - i <= 0)
                    {
                        return;
                    }
                    Matchresult *= (firstItem.size() - i);
                }
                if(type == pos){
                    num_positive_results_ += Matchresult;
                }
                else{
                    if(type == neg){
                        num_negative_results_ += Matchresult;
                    }
                    else{
                        num_initial_results_ += Matchresult;
                    }
                }
                
                return;
            }
        }
        std::vector<std::vector<int>> needToCombine;
        for(int i = 0; i < needToCombineV1.size(); i++){
            const auto & I_isolateVertexCandidate = needToCombineV1[i];
            std::vector<int> I_needToCombine;
            I_needToCombine.reserve(I_isolateVertexCandidate.size());
            int I_NoOverLeafWeight = 0;
            for(int k = 0; k < I_isolateVertexCandidate.size(); k++){
                int vertex = I_isolateVertexCandidate[k];
                if(isolatedVertexMap[vertex] > 1){
                    I_needToCombine.push_back(vertex);
                }
                else{
                    I_NoOverLeafWeight++;
                }
            }
            NoOverLeafWeight.push_back(I_NoOverLeafWeight);
            if(I_NoOverLeafWeight > 0){
                I_needToCombine.push_back(-1);
            }
            needToCombine.push_back(I_needToCombine);
        }
        DFS_SEARCH_ALL(needToCombine, NoOverLeafWeight, isolateVertexIndex, 0, type, 1);
    }
}

/**
 * @description: DFS match
 * @param {vector<uint>} &
 * @param {uint} idx
 * @param {searchType} type
 * @return {*}
 */
void CSMPP::DFS_SEARCH_ALL(const std::vector<uint> & needToCombine, const std::vector<uint> & cacheMap, uint idx, searchType type)
{
    uint posCurrent = needToCombine[idx];
    for(auto & v : this->matchCandidate[cacheMap[posCurrent]])
    {
        if(visited_[v]) continue;
        matchVertex(v, posCurrent);
        if(idx + 1 < needToCombine.size())
        {
            DFS_SEARCH_ALL(needToCombine, cacheMap, idx + 1, type);
        }
        else
        {
            this->printMatch();
            if(type == pos)
            {
                num_positive_results_++;
            }
            else
            {
                if(type == neg)
                {
                    num_negative_results_++;
                }
                else
                {
                    this->num_initial_results_++;
                }
            }
        }
        popVertex(v, posCurrent);
    }
}

/**
 * @description: DFS match use weight
 * @param {vector<std::vector<int>>} &
 * @param {vector<uint>} &
 * @param {vector<uint>} &
 * @param {uint} idx
 * @param {searchType} type
 * @param {size_t} weight
 * @return {*}
 */
void CSMPP::DFS_SEARCH_ALL(const std::vector<std::vector<int>> & needToCombine, const std::vector<uint> & NoOverLeafWeight, const std::vector<uint> & isolatedIndex, uint idx, searchType type, size_t weight)
{
    for(auto & v : needToCombine[idx])
    {
        if(v != -1 && visited_[v])continue;
        if(v == -1)
            weight *= NoOverLeafWeight[idx];
        else
            matchVertex(v, isolatedIndex[idx]);
        if(idx + 1 < needToCombine.size())
        {
            DFS_SEARCH_ALL(needToCombine, NoOverLeafWeight, isolatedIndex, idx + 1, type, weight);
        }
        else
        {
            if(type == pos)
            {
                num_positive_results_ += weight;
            }
            else
            {
                if(type == neg)
                {
                    num_negative_results_ += weight;
                }
                else
                {
                    this->num_initial_results_ += weight;
                }
            }
        }
        if(v == -1)
            weight /= NoOverLeafWeight[idx];
        else
            popVertex(v, isolatedIndex[idx]);
    }
}

/**
 * @description: print final result
 * @return {*}
 */
void CSMPP::printMatch()
{
    std::cout << "print match: " << std::endl;
    for(auto & u : this->match)
    {
        std::cout << u << " ";
    }
    std::cout << std::endl;
}
#endif