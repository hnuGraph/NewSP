#ifndef UTILS_TYPES
#define UTILS_TYPES

#include <chrono>
#include <iostream>
#include <climits>
#include <functional>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <stack>
#include <algorithm>   // for std::shuffle
#include <random>      // for std::default_random_engine
#include <iterator>    // for std::ostream_iterator
#include "assert.h"


#define NOT_EXIST INT_MAX
#define UNMATCHED INT_MAX

enum vertexType
{
    freeVertex, freezeVertex, isolatedVertex
};
enum LRAndIndexCheckType {
    Part1Check/*do just check without use lr in part 1*/,
    Part2Check/*do just check without use lr in part 2*/,
    Part1Nothing/*do nothing in part 1*/,
    Part2Nothing/*do nothing in part 2*/,
    ioslatedVertexNothing/*isolated vertex nothing to do*/
};
enum unFreezeStackType { 
    finalStack, runningStack
};
enum useCacheStatus
{
    useAndFinsh, useNotFinsh, notUse
};
class Timer {

public:

    Timer() : status(TimerStatus::off), time_cost_ms_(0) {}

    enum class TimerStatus {
        on, off
    };

    void StartTimer() {
        assert(status == TimerStatus::off);
        status = TimerStatus::on;
        time_start_ = std::chrono::steady_clock::now();

    };


    void StopTimer() {
        assert(status == TimerStatus::on);
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        status = TimerStatus::off;
        time_cost_ms_ += std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_start_);
    };

    double countUillisSecond() {

    }

    double CountTime() {
//        return (double)time_cost_ms_.count()*std::chrono::milliseconds::period::num/std::chrono::milliseconds::period::den;
        return (double)time_cost_ms_.count();
//        return time_cost_ms_.count();
    }
    void ResetTime(){
        time_cost_ms_=0*std::chrono::nanoseconds(1);
    }

    bool oneMinTest(){
        StopTimer();
        //std::cout << "time over time: " << (double)time_cost_ms_.count() << std::endl;
        if((double)time_cost_ms_.count() > 60000000000.00){
            std::cout << "time over time: " << (double)time_cost_ms_.count() << std::endl;
            return true;
        }
        StartTimer();
        return false;
    }

private:
    TimerStatus status;
    std::chrono::steady_clock::time_point time_start_;
    std::chrono::nanoseconds time_cost_ms_;
};

// Time counting
#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast<\
    std::chrono::nanoseconds>(Get_Time() - start).count()/(float)1000
#define Print_Time(str, start) std::cout << std::fixed << std::setprecision(7)<< str << Duration(start) << std::endl

struct InsertUnit {
    char type;  // 'v' or 'e' 
    bool is_add;// addition or deletion
    uint id1;   // vertex id or edge source id
    uint id2;   // edge target id
    uint label; // vertex or edge label
    InsertUnit(char type_arg, bool is_add_arg, uint id1_arg, uint id2_arg, uint label_arg)
    : type(type_arg), is_add(is_add_arg), id1(id1_arg), id2(id2_arg), label(label_arg) {}
};

class sampleUnit{
public:
    int minSize = 0;
    double avgDegreeSize = 0;
    int rank = 0;
    int motifSize = 0;
    int sizeMinusDepth = 0;
    double sizeMinusDepthDivisionSize = 0;
    int descNeighborMinus1 = 0;
    int freezeVertexAfter = 0;
    int isolatedVertexAfter = 0;
    vertexType currentVertexType = freeVertex;
    Timer motifTimer, nomotifTimer;
};

class fileSystem{
    private:
        std::string sampleSavePath;
        std::string trainFilePath;
        std::string predictionFilePath;
        std::string validationFilePath;
        std::ifstream readSampleFile;
        std::ifstream readTrainFile;
        std::ifstream readPredictionFile;
        std::ifstream readValidationFile;
        std::ofstream writeSampleFile;
        std::ofstream writeTrainFile;
        std::ofstream writePredictionFile;
        std::ofstream writeVaildationFile;
    public:
        void init(std::string sampleSavePath){
            this->sampleSavePath = sampleSavePath;
        }
        void writeSampleOpen(){
            this->writeSampleFile.open(this->sampleSavePath, std::ios::out);
        }
        void writeSampleToFile(std::string sampleStr){
            this->writeSampleFile << sampleStr;
        }
        void writeSampleClose(){
            this->writeSampleFile.close();
        }
};

class sample{
private:
    int Threshold;
    std::string filePath2SEE;
    std::string filePath;
    int times;
    bool recordBotton;
    std::ofstream fp;
    std::ofstream fp2SEE;
    uint PosSample;
    uint NegSample;

public:
    void setSample(int Treshold, std::string filePath, std::string filePath2See){
        this->Threshold = Treshold;
        this->filePath = filePath;
        this->filePath2SEE = filePath2See;
        this->times = 0;
        this->recordBotton = true;
        fp.open(this->filePath, std::ios::out);
        fp2SEE.open(this->filePath2SEE, std::ios::out);
        this->PosSample = 0;
        this->NegSample = 0;
    }
    ~sample(){
        fp.close();
        fp2SEE.close();
    }
    std::string sampleUnitToStr(sampleUnit& _sampleUnit){
        std::string sampleStr = "";
        double time1 = _sampleUnit.motifTimer.CountTime()/1000.0;
        double time0 = _sampleUnit.nomotifTimer.CountTime()/1000.0;
        double timeDiff = 0.0;
        if(time1 < time0){
            sampleStr += "1 ";
            timeDiff = time0 - time1;
            this->PosSample++;
        }
        else{
            sampleStr += "0 ";
            timeDiff = time1 - time0;
            this->NegSample++;
        }
        sampleStr += std::to_string(time1) + " ";
        sampleStr += std::to_string(time0) + " ";
        sampleStr += std::to_string(timeDiff) + " ";
        sampleStr += std::to_string(_sampleUnit.rank) + " ";
        sampleStr += std::to_string(_sampleUnit.sizeMinusDepth) + " ";
        sampleStr += std::to_string(_sampleUnit.sizeMinusDepthDivisionSize) + " ";
        sampleStr += std::to_string(_sampleUnit.descNeighborMinus1) + " ";
        sampleStr += std::to_string(_sampleUnit.freezeVertexAfter) + " ";
        sampleStr += std::to_string(_sampleUnit.isolatedVertexAfter) + " ";
        sampleStr += std::to_string(_sampleUnit.minSize) + " ";
        sampleStr += std::to_string(_sampleUnit.avgDegreeSize) + "\n";
        return sampleStr;
    }
    void setRecordButton(bool status){
        this->recordBotton = status;
    }
    bool getRecordButton(){
        return this->recordBotton;
    }
    uint getPosSampleCount(){
        return this->PosSample;
    }
    uint getNegSampleCount(){
        return this->NegSample;
    }
};

// from boost (functional/hash):
// see http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html template
template <typename T>
inline void hash_combine(std::size_t &seed, const T &val) {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
// auxiliary generic functions to create a hash value using a seed
template <typename T> inline void hash_val(std::size_t &seed, const T &val) {
    hash_combine(seed, val);
}
template <typename T, typename... Types>
inline void hash_val(std::size_t &seed, const T &val, const Types &... args) {
    hash_combine(seed, val);
    hash_val(seed, args...);
}

template <typename... Types>
inline std::size_t hash_val(const Types &... args) {
    std::size_t seed = 0;
    hash_val(seed, args...);
    return seed;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        return hash_val(p.first, p.second);
    }
};


#endif //UTILS_TYPES
