#ifndef GRAPHENE_HPP
#define GRAPHENE_HPP

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class Graphene{
    typedef matrix<complex_t> state_type;
public:
    virtual void propagate(
        const state_type& rho, state_type& drhodt, const double t,
        const double& kx_t, const double& ky_t) const=0;

    virtual size_t nstates() const=0;
private:
};

#endif