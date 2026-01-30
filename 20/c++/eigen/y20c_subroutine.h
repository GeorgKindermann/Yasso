#ifndef YASSO_H
#define YASSO_H

#include <array>
#include <memory>
#include <vector>

namespace yasso {

class yasso20 {
public:
    yasso20();
    yasso20(const std::array<double, 35> &theta);
    ~yasso20();

    void setTheta(const std::array<double, 35> &theta);
    void setClimSizeLeach(const std::array<double, 12>& avgT, double sumP, double diam, double leach);
    void setTimespan(double timespan);
    
    bool isThereDecomposition() const { return !noDecomposition; }
    
    void getSpin(const std::array<double, 5>& infall, std::array<double, 5>& result, int years = 700);
    void getNextTimestep(const std::array<double, 5>& init, const std::array<double, 5>& infall, std::array<double, 5>& result);
    
    size_t setTaylorTerms(size_t terms);

private:
    struct Impl; 
    std::unique_ptr<Impl> pimpl;

    double timespan = 1.0;
    bool noDecomposition = true;
    size_t taylorTerms = 11;
};

} // namespace yasso

#endif
