#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "graph.h"
/**
 * @description: add vertex for graph
 * @param {uint} id
 * @param {uint} label
 * @return {*}
 */
void Graph::AddVertex(uint id, uint label)
{
    if (id >= vlabels_.size())
    {
        vlabels_.resize(id + 1, NOT_EXIST);
        vlabels_[id] = label;
        neighbors_.resize(id + 1);
        elabels_.resize(id + 1);
    }
    else if (vlabels_[id] == NOT_EXIST)
    {
        vlabels_[id] = label;
    }
    
    vlabel_count_ = std::max(vlabel_count_, label + 1);
}

/**
 * @description: remove the query vertex
 * @param {uint} id
 * @return {*}
 */
void Graph::RemoveVertex(uint id)
{
    vlabels_[id] = NOT_EXIST;
    neighbors_[id].clear();
    elabels_[id].clear();
}

/**
 * @description: add edge for graph
 * @param {uint} v1
 * @param {uint} v2
 * @param {uint} label
 * @return {*}
 */
void Graph::AddEdge(uint v1, uint v2, uint label)
{

    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower != neighbors_[v1].end() && *lower == v2) return;
    size_t dis = std::distance(neighbors_[v1].begin(), lower);
    neighbors_[v1].insert(lower, v2);
    elabels_[v1].insert(elabels_[v1].begin() + dis, label);
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].insert(lower, v1);
    elabels_[v2].insert(elabels_[v2].begin() + dis,label);
    if(v1 > v2)
        std::swap(v1, v2);
    this->edge.emplace_back(Edge(v1, v2, this->GetVertexLabel(v1), this->GetVertexLabel(v2), label, this->edge.size()));
    edge_count_++;
    elabel_count_ = std::max(elabel_count_, label + 1);
}

/**
 * @description: remove graph's edge 
 * @param {uint} v1
 * @param {uint} v2
 * @return {*}
 */
void Graph::RemoveEdge(uint v1, uint v2)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower == neighbors_[v1].end() || *lower != v2)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v1].erase(lower);
    elabels_[v1].erase(elabels_[v1].begin() + std::distance(neighbors_[v1].begin(), lower));
    
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);
    if (lower == neighbors_[v2].end() || *lower != v1)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v2].erase(lower);
    elabels_[v2].erase(elabels_[v2].begin() + std::distance(neighbors_[v2].begin(), lower));
}

uint Graph::GetVertexLabel(uint u) const
{
    return vlabels_[u];
}

const std::vector<uint>& Graph::GetNeighbors(uint v) const
{
    return neighbors_[v];
}

const std::vector<uint>& Graph::GetNeighborLabels(uint v) const
{
    return elabels_[v];
}

std::tuple<uint, uint, uint> Graph::GetEdgeLabel(uint v1, uint v2) const
{
    uint v1_label, v2_label, e_label;
    v1_label = GetVertexLabel(v1);
    v2_label = GetVertexLabel(v2);

    const std::vector<uint> *nbrs;
    const std::vector<uint> *elabel;
    uint other;
    if (GetDegree(v1) < GetDegree(v2))
    {
        nbrs = &GetNeighbors(v1);
        elabel = &elabels_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        elabel = &elabels_[v2];
        other = v1;
    }
    
    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_label = elabel->at(mid);
            return {v1_label, v2_label, e_label};
        }
    }
    return {v1_label, v2_label, -1};
}

uint Graph::GetDegree(uint v) const
{
    return neighbors_[v].size();
}

uint Graph::GetDiameter() const
{
    uint diameter = 0;
    for (uint i = 0u; i < NumVertices(); i++)
    if (GetVertexLabel(i) != NOT_EXIST)
    {
        std::queue<uint> bfs_queue;
        std::vector<bool> visited(NumVertices(), false);
        uint level = UINT_MAX;
        bfs_queue.push(i);
        visited[i] = true;
        while (!bfs_queue.empty())
        {
            level++;
            uint size = bfs_queue.size();
            for (uint j = 0u; j < size; j++)
            {
                uint front = bfs_queue.front();
                bfs_queue.pop();

                const auto& nbrs = GetNeighbors(front);
                for (const uint nbr: nbrs)
                {
                    if (!visited[nbr])
                    {
                        bfs_queue.push(nbr);
                        visited[nbr] = true;
                    }
                }
            }
        }
        if (level > diameter) diameter = level;
    }
    return diameter;
}

/**
 * @description: load singel query graph
 * @param {string} &path
 * @return {*}
 */
void Graph::LoadFromFile(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    char type;
    while (ifs >> type)
    {
        if (type == 't')
        {
            char temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            AddVertex(vertex_id, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            AddEdge(from_id, to_id, label);
        }
    }
    ifs.close();
}

/**
 * @description: upload update stream
 * @param {string} &path
 * @return {*}
 */
void Graph::LoadUpdateStream(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    std::string type;
    while (ifs >> type)
    {
        if (type == "v" || type == "-v")
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            updates_.emplace('v', type == "v", vertex_id, 0u, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            updates_.emplace('e', type == "e", from_id, to_id, label);
        }
    }
    ifs.close();
    std::cout << "# update size = " << this->updates_.size() << std::endl;
}

void Graph::PrintMetaData() const
{
    std::cout << "# vertices = " << NumVertices() << 
        "  # edges = " << NumEdges() << 
        "  vlabel size = " << this->NumVLabels() << 
        "  elabel size = " << this->NumELabels() << std::endl;
}

const std::vector<uint> & Graph::GetMatchOrder(uint index) const{
    return this->matchOrder[index];
}
/**
 * @description: begin build matching order
 * @return {*}
 */
void Graph::MatchOrderInit(){
    uint numVertices = this->NumVertices();
    for(uint i = 0; i < this->edge.size(); i++){
        auto edge = this->edge[i];
        std::vector<bool>visit(numVertices, false);
        std::vector<uint> order{};
        EdgeMakeOrder(i, edge.GetV1(), edge.GetV2(), visit, order);
        this->matchOrder.push_back(order);
    }
}

/**
 * @description: make dependent list(left label)
 * @param {uint} edgeIndex
 * @param {vector<uint>} &
 * @param {uint} visitedVertex
 * @return {*}
 */
void Graph::DescListInit(uint edgeIndex, std::vector<uint> & order, uint visitedVertex){
    std::vector<std::tuple<uint,uint,uint>> neighborBefore{};
    const std::vector<uint> & VertexNeighbor = this->GetNeighbors(visitedVertex);
    for(int i = 0; i < order.size(); i++){
        if(find(VertexNeighbor.begin(), VertexNeighbor.end(), order[i]) != VertexNeighbor.end()){
            uint vertexLabel = this->GetVertexLabel(order[i]);
            uint eLabel = std::get<2>(this->GetEdgeLabel(visitedVertex, order[i]));
            neighborBefore.push_back(std::tuple(i, vertexLabel, eLabel));
        }
    }
    if(this->descList.size() == edgeIndex){
        std::vector<std::vector<std::tuple<uint,uint,uint>>> edgeIndexMatchDesList {};
        edgeIndexMatchDesList.push_back(neighborBefore);
        this->descList.push_back(edgeIndexMatchDesList);
    }
    else{
        this->descList[edgeIndex].push_back(neighborBefore);
    }
}

/**
 * @description: make matching order
 * @param {uint} edgeIndex
 * @param {uint} v1
 * @param {uint} v2
 * @param {vector<bool>} &
 * @param {vector<uint>} &
 * @return {*}
 */
void Graph::EdgeMakeOrder(uint edgeIndex, uint v1, uint v2, std::vector<bool> & visit, std::vector<uint> & order){
    uint vMax = this->NumVertices();
    if(order.size() == vMax){
        return;
    } 
    if(order.size() < 2){
        if(this->GetDegree(v1) > this->GetDegree(v2)){
            std::swap(v1, v2);
        }
        this->DescListInit(edgeIndex, order, v1);
        order.push_back(v1);
        visit[v1] = true;

        this->DescListInit(edgeIndex, order, v2);
        order.push_back(v2);
        visit[v2] = true;
        
        EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
    }
    else{
        //round 1 check neighbor has been visited
        std::map<uint,uint>compareUse;
        for(uint i = 0; i < vMax; i++){
            if(visit[i] == false){
                std::vector<uint>neighbor = this->GetNeighbors(i);
                for(auto item : order){
                    if(find(neighbor.begin(), neighbor.end(), item) != neighbor.end()){
                        compareUse[i]++;
                    }
                }
            }
        }
        std::vector<std::pair<uint,uint> > arr;
        for (std::map<uint, uint>::iterator it=compareUse.begin();it!=compareUse.end();++it)
        {
            arr.push_back(std::make_pair(it->first,it->second));
        }
        sort(arr.begin(),arr.end(),cmp);
        std::vector<uint>Candidate;
        uint max = 0;
        for(auto item : arr){
            if(item.second >= max){
                max = item.second;
                Candidate.push_back(item.first);
            }
            else break;
        }
        if(Candidate.size() == 1){
            this->DescListInit(edgeIndex, order, Candidate[0]);
            visit[Candidate[0]] = true;
            order.push_back(Candidate[0]);
            EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
            return;
        }
        //round 2 check degree
        compareUse.clear();
        arr.clear();
        max = 0;
        for(uint i = 0; i < Candidate.size(); i++){
            compareUse[Candidate[i]] = this->GetDegree(Candidate[i]);
        }
        for (std::map<uint, uint>::iterator it=compareUse.begin();it!=compareUse.end();++it)
        {
            arr.push_back(std::make_pair(it->first,it->second));
        }
        sort(arr.begin(),arr.end(),cmp);
        Candidate.clear();
        max = 0;
        for(auto item : arr){
            if(item.second >= max){
                max = item.second;
                Candidate.push_back(item.first);
            }
            else break;
        }
        if(Candidate.size() == 1){
            this->DescListInit(edgeIndex, order, Candidate[0]);
            visit[Candidate[0]] = true;
            order.push_back(Candidate[0]);
            EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
            return;
        }
        //round 3 check id
        sort(Candidate.begin(), Candidate.end());
        this->DescListInit(edgeIndex, order, Candidate[0]);
        visit[Candidate[0]] = true;
        order.push_back(Candidate[0]);
        EdgeMakeOrder(edgeIndex, v1, v2, visit, order);
        return;
    }
    
}

/**
 * @description: get matching order by edge
 * @param {Edge} edge
 * @return {*}
 */
const std::vector<uint> & Graph::EdgeGetOrder(Edge edge) const{
    return this->matchOrder[edge.GetIndex()];
}

/**
 * @description: get matching order by edge
 * @param {uint} index
 * @return {*}
 */
const std::vector<uint> & Graph::EdgeGetOrder(uint index) const{
    return this->matchOrder[index];
}

//build query indenpendent index
void Graph::indexInit(){
#if GRAPH_TYPE == 0
    uint vLabelMax = this->NumVLabels();
    uint eLabelMax = this->NumELabels();
    uint vMax = this->NumVertices();
    this->index.resize(vMax, nullptr);
    for(uint k = 0; k < vMax; k++){
        std::vector<uint>iNeighbor = this->GetNeighbors(k);
        std::vector<uint>iNeighborEdgeLabel = this->GetNeighborLabels(k);
        int * indexDistribution = new int [vLabelMax + eLabelMax]();
        for(int i = 0; i < iNeighbor.size(); i++){
            indexDistribution[this->GetVertexLabel(iNeighbor[i])]++;
            indexDistribution[vLabelMax + iNeighborEdgeLabel[i]]++;
        }
        this->index[k] = indexDistribution;
    }
    //for dense graph
#elif GRAPH_TYPE == 1
    this->index.resize(this->NumVertices());
    for(uint VertexIndex = 0; VertexIndex < this->NumVertices(); VertexIndex++)
    {
        const std::vector<uint> & iNeighbor = this->GetNeighbors(VertexIndex);
        const std::vector<uint> & iNeighborEdgeLabel = this->GetNeighborLabels(VertexIndex);
        for(auto item : iNeighbor)
        {
            this->index[VertexIndex][this->GetVertexLabel(item)]++;
        }
        for(auto item : iNeighborEdgeLabel)
        {
            this->index[VertexIndex][this->NumVLabels() + item]++;
        }
    }
#endif
}

/**
 * @description: update query-independent index
 * @param {uint} v1
 * @param {uint} v2
 * @param {uint} label
 * @param {bool} op
 * @return {*}
 */
void Graph::indexUpdate(uint v1, uint v2, uint label, bool op){
    uint v1Label = this->GetVertexLabel(v1);
    uint v2Label = this->GetVertexLabel(v2);
    uint Label = this->NumVLabels() + label;
    if(op == true){
        this->index[v1][v2Label]++;
        this->index[v1][Label]++;
        this->index[v2][v1Label]++;
        this->index[v2][Label]++;
    }
    else{
        this->index[v1][v2Label]--;
        this->index[v1][Label]--;
        this->index[v2][v1Label]--;
        this->index[v2][Label]--;
    }
}

Edge Graph::GetEdge(uint index){
    return this->edge[index];
}

std::vector<Edge> Graph::GetEdge(){
    return this->edge;
}

/**
 * @description: mapping vertex
 * @param {uint} vertex
 * @param {uint} beginEdge
 * @return {*}
 */
bool Graph::mappingAdd(uint vertex, uint beginEdge){
    this->mapping.push_back(vertex);
    std::cout << "# Mapping Add size is : " << this->mapping.size() << std::endl;
    if(this->mapping.size() == this->NumVertices()){
        std::cout << "match success" << std::endl;
        for(uint i = 0; i < this->NumVertices(); i++){
            std::cout << "query vertex: " << this->GetMatchOrder(beginEdge)[i] << "--> data graph: " << this->mapping[i] << std::endl;
        }
        return true;
    }
    return false;
}

void Graph::mappingPopTail(){
    this->mapping.pop_back();
}

//desclist
std::vector<std::tuple<uint, uint, uint>> Graph::GetDescList(uint edgeIndex, uint depth){
    return this->descList[edgeIndex][depth];
}

//load batch query graph
void Graph::LoadFromFile(const std::string &prefixPath, std::vector<Graph> & queryGraph){
    std::vector<std::string>fileNames;
    struct stat s;
    stat(prefixPath.c_str(), &s);
    DIR* open_dir = opendir(prefixPath.c_str());
    if (NULL == open_dir) {
        std::exit(EXIT_FAILURE);
    }
    dirent* p = nullptr;
    while( (p = readdir(open_dir)) != nullptr) {
        struct stat st;
        if (p->d_name[0] != '.') {
            std::string name = prefixPath + std::string("/") + std::string(p->d_name);
            stat(name.c_str(), &st);
            fileNames.push_back(name);
        }
    }
    closedir(open_dir);
    for(int i = 0; i < fileNames.size(); i++){
        Graph queryGraphItem;
        queryGraphItem.LoadFromFile(fileNames[i]);
        queryGraph.push_back(queryGraphItem);
    }
    
    std::cout << "Load query graph number: " << queryGraph.size() << std::endl;
}
//for test
void Graph::indexPrint(){
    int vLabelMax = this->NumVLabels();
    int eLabelMax = this->NumELabels();
    int vMAx = this->NumVertices();
    for(int i = 0; i < vMAx; i++){
        std::cout << "vLabel: ";
        for(int k = 0; k < vLabelMax; k++){
            std::cout << this->index[i][k] << " ";
        }
        std::cout << std::endl;
        std::cout << "eLabel: ";
        for(int k = 0; k < eLabelMax; k++){
            std::cout << this->index[i][k + vLabelMax] << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;

    }
}
//for test
void Graph::MatchOrderAllPrint(){
    for(int i = 0; i < this->NumEdges(); i++){
        std::cout << "egde index: " << i << std::endl;
        // i is edge index
        const auto & matchOrder = this->matchOrder[i];
        const auto & vertexTypes = this->matchVertexTypes[i];
        const auto & desitemList = this->descList[i];
        const auto & unfreezeRecord_ = this->unfreezeRecord[i];
        const auto & cacheDependent = this->matchCacheOrder[i];
        const auto & cacheSaveBool = this->matchCacheSave[i];
        //2.cout a match order
        for(int k = 0; k < matchOrder.size(); k++){
            std::cout << matchOrder[k] << " ";
        }
        std::cout << std::endl;
        for(int k = 0; k < matchOrder.size(); k++){
            std::cout << "vertex : " << matchOrder[k] << std::endl;
            std::cout << "vertex index: " << k << std::endl;
            std::cout << "vertex type: " << vertexTypes[k] << std::endl;
            std::cout << "cache dependent: " << cacheDependent[k] << std::endl;
            std::cout << "cache save: " << cacheSaveBool[k] << std::endl;
            const auto & unfreezerecord_ = unfreezeRecord_[k];
            std::cout << "unfreeze list : ";
            for(auto vertex : unfreezerecord_){
                std::cout << "[" << vertex.first << ", " << vertex.second << "]";
            }
            std::cout << std::endl;
            const auto & desitemList_ = desitemList[k];
            std::cout << "desItem list : " << std::endl;
            for(auto item : desitemList_){
                std::cout << "index: " << std::get<0>(item) << " item order index: " << this->matchOrder[i][std::get<0>(item)] << " item label: " << std::get<1>(item) << " elabel: " << std::get<2>(item) << std::endl ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

//for test
void Graph::DescListPrint(uint edgeIndex){
    std::vector<std::vector<std::tuple<uint,uint,uint>>> temp = this->descList[edgeIndex];
    std::cout << "desList size = " << temp.size() << std::endl;
    auto matchOrder = this->GetMatchOrder(edgeIndex);
    for(int index = 0; index < this->NumVertices(); index++){
        std::vector<std::tuple<uint,uint,uint>> temp_ = temp[index];
        std::cout << "match order index: " << index << " item is " << matchOrder[index] << " descList size is " << temp_.size() << std::endl;
        for(auto item : temp_){
            std::cout << "item in order : " << this->matchOrder[edgeIndex][std::get<0>(item)] << " item label:" << std::get<1>(item) << " elabel: " << std::get<2>(item) << std::endl ;
        }
        std::cout << std::endl;
    }
}

uint Graph::GetMappingSize(){
    return this->mapping.size();
}

void Graph::PrintGraphNature() const{
    double avgDegree = 0;
    uint maxDegree = 0;
    uint DegreeBiggerAvg = 0;
    const uint vMax = this->NumVertices();
    for(int i = 0; i < vMax; i++){
        uint tempDegree = this->GetDegree(i);
        avgDegree  += tempDegree;
        maxDegree = maxDegree > tempDegree? maxDegree : tempDegree;
    }
    avgDegree = avgDegree / vMax;
    for(int i = 0; i < vMax; i++){
        if(this->GetDegree(i) > avgDegree){
            DegreeBiggerAvg++;
        }
    }
    std::cout << "# avgDegree = " << avgDegree << std::endl;
    std::cout << "# maxDegree = " << maxDegree << std::endl;
    std::cout << "# DegreeBiggerAvg++ = " << DegreeBiggerAvg << std::endl;
}

void Graph::indexAvg() const{
#if GRAPH_TYPE == 0
    uint vLabelMax = this->NumVLabels();
    uint eLabelMax = this->NumELabels();
    const uint vMax = this->NumVertices();
    double * avg = new double [vLabelMax + eLabelMax]();
    for(int i =  0; i < vMax; i++){
        for(int k = 0; k < vLabelMax + eLabelMax; k++){
            avg[k] = avg[k] + this->index[i][k] ;
        }
    }
    std::cout << "avg index dirtribution: ";
    for(int k = 0; k < vLabelMax + eLabelMax; k++){
        avg[k] = avg[k] / vMax;
        std::cout << avg[k] << " ";
    }
    std::cout << std::endl;
#endif
}

const uint Graph::getIndexValue(uint v, uint pos) const{
#if GRAPH_TYPE == 0
    return this->index[v][pos];
#elif GRAPH_TYPE == 1
    auto iter = this->index[v].find(pos);
    if(iter != this->index[v].end())
    {
        return iter->second;
    }
    return UINT_MAX;
#endif
}

/**
 * @description: make vertex attri
 * @return {*}
 */
void Graph::matchOrderTypeSet(){
    // set vertexTypes
    for(int edgeIndex = 0; edgeIndex < this->NumEdges(); edgeIndex++){
        const auto & matchOrder = this->matchOrder[edgeIndex];
        const auto & desitemList = this->descList[edgeIndex];
        //0. set unfreeze record;
        std::vector<std::vector<std::pair<uint,uint>>> currentMatchOrderUnfreezeRecord;
        currentMatchOrderUnfreezeRecord.resize(matchOrder.size(), {});
        //1. get vertexTypes and get unfreeze record
        std::vector<vertexType> currentMatchOrderTypes;
        for(int i = 0; i < matchOrder.size(); i++){
            vertexType iType = isolatedVertex;
            if(i == 0 || i == 1){
                iType = freeVertex;
            }
            else{
                for(int k = i + 1; k < matchOrder.size(); k++){
                    const auto & desitemListTuples = desitemList[k];
                    for(int j = 0; j < desitemListTuples.size(); j++){
                        if(std::get<0>(desitemListTuples[j]) == i){
                            if(k == i + 1){
                                iType = freeVertex;
                            }
                            if(k != i + 1){
                                currentMatchOrderUnfreezeRecord[k].push_back(std::make_pair(i, std::get<2>(desitemListTuples[j]))); // not the after vertex of i , should record it can unfreeze which vertices
                                iType = freezeVertex;
                            }
                            break;
                        }
                    }
                    if(iType == freeVertex || iType == freezeVertex){
                        break;
                    }
                }
            }
            currentMatchOrderTypes.push_back(iType);
        }
        this->matchVertexTypes.push_back(currentMatchOrderTypes);
        this->unfreezeRecord.push_back(currentMatchOrderUnfreezeRecord);
    }
    //delete freeze vertex from unfreeze vertex in desItemList
    for(int edgeIndex = 0; edgeIndex < this->NumEdges(); edgeIndex++){// for every match order
        const auto & unfreezeRecord_ = this->unfreezeRecord[edgeIndex];// get match order's unfreezeRecord
        for(int i = 0; i < unfreezeRecord_.size(); i++){
            if(unfreezeRecord_[i].size() != 0){//now unfreezeRecord_[i] is a unfreeze vertex
                auto & desitemItemOfi = this->descList[edgeIndex][i];// get unfreezeRecord_[i]'s desitem
                for(int k = 0; k < unfreezeRecord_[i].size(); k++){// loop of unfreezeRecord_[i]'s unfreeze record
                    auto iter = desitemItemOfi.begin();
                    while(iter != desitemItemOfi.end()){
                        if(std::get<0>(*iter) == unfreezeRecord_[i][k].first){
                            desitemItemOfi.erase(iter);//delete freeze vertex from unfreeze vertex 's desitemList
                            break;
                        }
                        iter++;
                    }
                }
            }
        }
    }
    //record isolate vertex
    for(int i = 0; i < this->NumEdges(); i++){
        const auto & vertexTypes = this->matchVertexTypes[i];
        std::vector<uint>isolateRecord;
        for(int k = 0; k < vertexTypes.size(); k++){
            if(vertexTypes[k] == isolatedVertex){
                isolateRecord.push_back(k);
            }
        }
        this->matchIsolatedOrder.push_back(isolateRecord);
    }
    //record freezeVertex + isolatedVertex Num After
    for(int i = 0; i < this->NumEdges(); i++){
        const auto & vertexTypes = this->matchVertexTypes[i];
        std::vector<uint> freezeVertexRecorder;
        std::vector<uint> isolatedVertexRecorder;
        freezeVertexRecorder.resize(vertexTypes.size());
        isolatedVertexRecorder.resize(vertexTypes.size());
        uint freezeCount = 0;
        uint isolateCount = 0;
        for(int k = vertexTypes.size() - 1; k >= 0; k--){
            freezeVertexRecorder[k] = freezeCount;
            isolatedVertexRecorder[k] = isolateCount;
            if(vertexTypes[k] == freezeVertex){
                freezeCount++;
            }
            if(vertexTypes[k] == isolatedVertex){
                isolateCount++;
            }
        }
        this->freezeVertexNumAfter.push_back(freezeVertexRecorder);
        this->isolatedVertexNumAfter.push_back(isolatedVertexRecorder);
    }
}

const vertexType Graph::getVertexType(uint edgeIndex, uint pos) const{
    return this->matchVertexTypes[edgeIndex][pos];
}

const std::vector<std::pair<uint, uint>> & Graph::getUnfreezeList(uint edgeIndex, uint pos)const {
    return this->unfreezeRecord[edgeIndex][pos];
}

void Graph::setVertexFree(uint edgeIndex, uint pos){
    this->matchVertexTypes[edgeIndex][pos] = freeVertex;
}

void Graph::setVertexStatus(uint edgeIndex, const std::vector<std::pair<uint,uint>> & vertexs, vertexType type){
    for(int i = 0; i < vertexs.size(); i++){
        this->matchVertexTypes[edgeIndex][vertexs[i].first] = type;
    }
}

const std::vector<uint> & Graph::getIsolatedVertexIndex(uint edgeIndex) const{
   return this->matchIsolatedOrder[edgeIndex];
}

const std::vector<vertexType> & Graph::getVertexType(uint edgeIndex) const{
    return this->matchVertexTypes[edgeIndex];
}

void Graph::isolatedVertexTimesAdd(const std::vector<uint> & candidateVertexs){
    for(int i = 0; i < candidateVertexs.size(); i++){
        this->isolatedVertexTimes[candidateVertexs[i]]++;
    }
}

void Graph::isolatedVertexTimesMinus(const std::vector<uint> & candidateVertexs){
    for(int i = 0; i < candidateVertexs.size(); i++){
        this->isolatedVertexTimes[candidateVertexs[i]]--;
    }
}

const uint Graph::getFreezeVertexNumAfter(uint edgeIndex, uint pos) const{
    return this->freezeVertexNumAfter[edgeIndex][pos];
}

const uint Graph::getIsolatedVertexNumAfter(uint edgeIndex, uint pos) const{
    return this->isolatedVertexNumAfter[edgeIndex][pos];
}

const uint Graph::getDestListSize(uint egdeIndex, uint depth) const{
    return this->descList[egdeIndex][depth].size();
}

const uint Graph::getUnfreezeListSize(uint edgeIndex, uint depth) const{
    return this->unfreezeRecord[edgeIndex][depth].size();
}

/**
 * @description: make cache use index
 * @return {*}
 */
void Graph::initCacheOrder()
{
    this->matchCacheOrder.resize(this->NumEdges());
    this->matchCacheSave.resize(this->NumEdges());
    this->cacheMap.resize(this->NumEdges());
    const uint LabelNum = this->NumELabels() + this->NumVLabels();
    for(int i = 0; i < this->NumEdges(); i++)
    {
        const auto & mOrder = this->GetMatchOrder(i);
        std::vector<uint> cacheUse;
        std::vector<uint> Map;
        for(uint i = 0; i < this->NumVertices(); i++)
        {
            cacheUse.emplace_back(i);
            Map.emplace_back(i);
        }
        std::vector<bool> savePending(this->NumVertices(), false);
        for(int currentNode = this->NumVertices() - 1; currentNode > 1; currentNode--)
        {

            const auto & _CDependentList = this->GetDescList(i, currentNode);
            std::map<uint, uint> elabelCheck;
            for(auto & item : _CDependentList)
            {
                elabelCheck[std::get<0>(item)] = std::get<2>(item);
            }
            auto & _ClabelDistribution = this->index[mOrder[currentNode]];
            for(int Node = 2; Node < currentNode; Node++)
            {
                bool same = true;
                vertexType nodeStatus = this->getVertexType(i, Node);
                if(nodeStatus == freeVertex) continue;
                if(
                    this->GetVertexLabel(mOrder[Node]) == this->GetVertexLabel(mOrder[currentNode]) 
                    // && elabelCheck.find(Node) != elabelCheck.end()
                )
                {
                    auto & _NlabelDistribution = this->index[mOrder[Node]];
                    #if GRAPH_TYPE == 0
                    for(int k = 0; k < LabelNum; k++)
                    {
                        if(_NlabelDistribution[k] > _ClabelDistribution[k])
                        {
                            same = false;
                            break;
                        }
                    }
                    #elif GRAPH_TYPE == 1
                    auto iter = _NlabelDistribution.begin();
                    while(iter != _NlabelDistribution.end())
                    {
                        if(_ClabelDistribution.find(iter->first) == _ClabelDistribution.end() || iter->second > _ClabelDistribution[iter->first])
                        {
                            same = false;
                            break;
                        }
                        iter++;
                    }
                    #endif
                    if(same)
                    {
                        const auto & _NDependentList = this->GetDescList(i, Node);
                        for(auto & item : _NDependentList)
                        {
                            auto iter = elabelCheck.find(std::get<0>(item));
                            if(iter == elabelCheck.end() || iter != elabelCheck.end() && iter->second != std::get<2>(item))
                            {
                                same = false;
                                break;
                            }
                        }
                        const auto & _NUnfreezeList = this->getUnfreezeList(i, Node);
                        for(auto & item : _NUnfreezeList)
                        {
                            auto iter = elabelCheck.find(std::get<0>(item));
                            if(iter == elabelCheck.end() || iter != elabelCheck.end() && iter->second != std::get<1>(item))
                            {
                                same = false;
                                break;
                            }
                        }
                    }
                }
                else
                {
                    same = false;
                }
                if(same)
                {
                    const auto & _NDependentList = this->GetDescList(i, Node);
                    const auto & _NUnfreezeList = this->getUnfreezeList(i, Node);
                    auto & depentlist = this->descList[i][currentNode];
                    auto & unfreezelist = this->unfreezeRecord[i][currentNode];
                    //case 2 test all
                    cacheUse[currentNode] = Node;
                    savePending[Node] = true;
                    for(auto item : _NDependentList)
                    {
                        uint index = std::get<0>(item);
                        for(int k = 0 ; k < depentlist.size(); k++)
                        {
                            if(std::get<0>(depentlist[k]) == index)
                            {
                                depentlist.erase(depentlist.begin() + k);
                                break;
                            }
                        }
                    }
                    for(auto item : _NUnfreezeList)
                    {
                        uint index = std::get<0>(item);
                        for(int k = 0 ; k < depentlist.size(); k++)
                        {
                            if(std::get<0>(depentlist[k]) == index)
                            {
                                depentlist.erase(depentlist.begin() + k);
                                break;
                            }
                        }
                    }
                    if(depentlist.size() == 0 && unfreezelist.size() == 0)
                        Map[currentNode] = Node;
                }
                break;
            }
        }
        this->matchCacheOrder[i] = std::move(cacheUse);
        this->matchCacheSave[i] = std::move(savePending);
        this->cacheMap[i] = std::move(Map);
    }
    
}

/**
 * @description: get current node use which cache
 * @param {uint} edgeIndex
 * @param {uint} depth
 * @return {*}
 */
uint Graph::getCacheStatu(uint edgeIndex, uint depth)
{
    return this->matchCacheOrder[edgeIndex][depth];
}

/**
 * @description: get current node whether should save
 * @param {uint} edgeIndex
 * @param {uint} depth
 * @return {*}
 */
bool Graph::getSaveStatu(uint edgeIndex, uint depth)
{
    return this->matchCacheSave[edgeIndex][depth];
}

/**
 * @description: get cache 
 * @param {uint} edgeIndex
 * @return {*}
 */
const std::vector<uint> & Graph::getCacheStatu(uint edgeIndex)
{
    return this->matchCacheOrder[edgeIndex];
}

/**
 * @description: adaptive check decision index make 
 * @return {*}
 */
void Graph::initDecision()
{
    this->DecisionList.resize(this->NumEdges());
    for(int EdgeIndex = 0; EdgeIndex < this->NumEdges(); EdgeIndex++)
    {
        std::vector<LRAndIndexCheckType> temp(this->NumVertices());
        for(int VertexIndex = 0; VertexIndex < this->NumVertices(); VertexIndex++)
        {
            const auto & VertexDesListSize = this->getDestListSize(EdgeIndex, VertexIndex);
            const auto & VertexFreezeListSize = this->getUnfreezeListSize(EdgeIndex, VertexIndex);
            auto VertexType = this->getVertexType(EdgeIndex, VertexIndex);
            temp[VertexIndex] = this->decisionMakeSystem.makeDecision(VertexType, VertexDesListSize, VertexFreezeListSize);
        }
        this->DecisionList[EdgeIndex] = std::move(temp);
    }
}

/**
 * @description: get adaptive decision
 * @param {uint} EdgeIndex
 * @param {uint} VertexIndex
 * @return {*}
 */
LRAndIndexCheckType Graph::getDecision(uint EdgeIndex, uint VertexIndex)
{
    return this->DecisionList[EdgeIndex][VertexIndex];
}

/**
 * @description: get cache map
 * @param {uint} egdeIndex
 * @param {uint} VertexIndex
 * @return {*}
 */
uint Graph::getCacheMap(uint egdeIndex, uint VertexIndex)
{
    return this->cacheMap[egdeIndex][VertexIndex];
}

const std::vector<uint> & Graph::getCacheMap(uint edgeIndex) const
{
    return this->cacheMap[edgeIndex];
}