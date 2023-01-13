#ifndef GRAPHENE_HPP
#define GRAPHENE_HPP

#include "utils/utils.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class Graphene{
protected:
    typedef matrix<complex_t> state_type;
    typedef vector<complex_t> vector_t;
public:
    virtual void propagate(
        const state_type& rho, state_type& drhodt, const double t,
        const double& kx_t, const double& ky_t) const=0;

    virtual size_t nstates() const=0;

    virtual double get_energy(const double& kx, const double& ky, const size_t& ist) const=0;
    virtual double get_dipole(const double& kx, const double& ky, const size_t& ist, const size_t& jst, const size_t& dir) const=0;
    virtual double get_energy_grad(const double& kx, const double& ky, const size_t& ist, const size_t& dir) const=0;

    virtual vector_t get_state(const double& kx, const double& ky, const size_t& ist) const=0;
private:
};

#endif