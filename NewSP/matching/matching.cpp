#include <unordered_set>
#include <vector>

#include <chrono>
#include <climits>
#include <functional>
#include <stdlib.h>
#include <iomanip>

#include "utils/types.h"
#include "graph/graph.h"
#include "matching/matching.h"

matching::matching(std::vector<Graph>&queryVec, Graph& query_graph, Graph& data_graph,
        size_t max_num_results, 
        bool print_prep,
        bool print_enum, 
        bool homo)
: query_(query_graph)
, data_(data_graph)
, queryVec(queryVec)

, max_num_results_(max_num_results)
, print_preprocessing_results_(print_prep)
, print_enumeration_results_(print_enum)
, homomorphism_(homo)

, visited_(data_.NumVertices(), false)
, num_initial_results_(0ul)
, num_positive_results_(0ul)
, num_negative_results_(0ul)
, num_intermediate_results_before_index_check_(0ul)
, num_intermediate_results_after_index_check_(0ul)
, num_intermediate_results_after_joinability_check_(0ul)
, num_intermediate_results_after_visit_check_(0ul)
, num_intermediate_results_with_empty_candidate_set_(0ul)
, num_intermediate_results_without_results_(0ul)
{}

void matching::Preprocessing()
{}

void matching::InitialMatching()
{}

void matching::AddEdge(uint v1, uint v2, uint label)
{
    data_.AddEdge(v1, v2, label);
}

void matching::RemoveEdge(uint v1, uint v2)
{
    data_.RemoveEdge(v1, v2);
}

void matching::AddVertex(uint id, uint label)
{
    data_.AddVertex(id, label);
}

void matching::RemoveVertex(uint id)
{
    data_.RemoveVertex(id);
}

void matching::GetMemoryCost(size_t &num_edges, size_t &num_vertices)
{
    num_edges = 0ul;
    num_vertices = 0ul;
}

void matching::GetNumInitialResults(size_t &num_initial_results)
{
    num_initial_results = num_initial_results_;
}

void matching::GetNumPositiveResults(size_t &num_positive_results)
{
    num_positive_results = num_positive_results_;
}

void matching::GetNumNegativeResults(size_t &num_negative_results)
{
    num_negative_results = num_negative_results_;
}

void matching::PrintCounter()
{
    std::cout << num_intermediate_results_before_index_check_ << " intermediate results before index check.\n";
    std::cout << num_intermediate_results_after_index_check_ << " intermediate results after index check.\n";
    std::cout << num_intermediate_results_after_joinability_check_ << " intermediate results after joinability check.\n";
    std::cout << num_intermediate_results_after_visit_check_ << " intermediate results after visit check.\n";
    std::cout << num_intermediate_results_with_empty_candidate_set_ << " intermediate results with empty candidate set.\n";
    std::cout << num_intermediate_results_without_results_ << " intermediate results without results.\n";
}

void matching::TimePrint(bool motif)
{
    /*
    std::cout << "Time Check" << std::endl;
    std::cout << "μs version (1μs = 1/1000000s)" << std::endl;
    std::cout << "# indexBuild time = " << (double)this->indexBuildTime * std::chrono::microseconds::period::num << "μs" << std::endl;
    std::cout << "# indexupdate time = " << (double)this->indexupdateTime * std::chrono::microseconds::period::num << "μs" << std::endl;
    std::cout << "# indexCheck time = " << (double)this->indexCheckTime * std::chrono::microseconds::period::num << "μs" << std::endl;

    std::cout << "# matchOrderBuild time = " << (double)this->matchOrderBuildTime * std::chrono::microseconds::period::num << "μs" << std::endl;
    std::cout << "# findQueryGraph time = " << (double)this->findQueryGraphTime * std::chrono::microseconds::period::num << "μs" << std::endl;

    std::cout << "# searchVertex time = " << (double)this->searchVertexTime * std::chrono::microseconds::period::num << "μs" << std::endl;

    std::cout << std::endl;
    std::cout << "s version" << std::endl;s
    std::cout << "# indexBuild time = " << (double)this->indexBuildTime * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << "s" << std::endl;
    std::cout << "# indexupdate time = " << (double)this->indexupdateTime * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << "s" << std::endl;
    std::cout << "# indexCheck time = " << (double)this->indexCheckTime * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << "s" << std::endl;

    std::cout << "# matchOrderBuild time = " << (double)this->matchOrderBuildTime * std::chrono::microseconds::period::num /std::chrono::microseconds::period::den << "s" << std::endl;
    std::cout << "# findQueryGraph time = " << (double)this->findQueryGraphTime * std::chrono::microseconds::period::num /std::chrono::microseconds::period::den<< "s" << std::endl;

    std::cout << "# searchVertex time = " << (double)this->searchVertexTime * std::chrono::microseconds::period::num/std::chrono::microseconds::period::den << "s" << std::endl;
    */
    //std::cout << "# desTime = " << (this->DescListTime)/(float)1000 << "ms" << stsd::endl;
    //std::cout << "# wcTime = " << (this->WCTime)/(float)1000 << "ms" << std::endl;
    if(motif){
        std::cout << std::fixed << std::setprecision(7)  << "# indexupdate = " << (this->indexupdateTime.CountTime())/(float)1000 << " microsecond" << std::endl;
        std::cout << std::fixed << std::setprecision(7)  << "# indexCheck = " << (this->indexCheckTime.CountTime())/(float)1000 << " microsecond" << std::endl;
        std::cout << std::fixed << std::setprecision(7)  << "# searchInitTime = " << (this->searchInitTime.CountTime())/(float)1000 << " microsecond" << std::endl;
    }
    std::cout << std::fixed << std::setprecision(7)  << "# WCTime = " << (this->WCTime.CountTime())/(float)1000 << " microsecond" << std::endl;
    std::cout << std::fixed << std::setprecision(7)  << "#LR predict = " << this->LRTime.CountTime()/1000.0 << " microsecond" << std::endl;
}