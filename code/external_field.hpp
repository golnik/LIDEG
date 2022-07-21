#ifndef EXTERNALFIELD_HPP
#define EXTERNALFIELD_HPP

class ExternalField{
public:
    virtual double operator()(const double& t) const=0;
private:
};

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using namespace boost::math::interpolators;

class ExternalFieldFromData:
public ExternalField{
    typedef cardinal_cubic_b_spline<double> spline_t;
public:
    ExternalFieldFromData(const std::vector<double>& data,
        const double& t0, const double& dt){
        _spl=new spline_t(data.begin(),data.end(),t0,dt);
    }
    ~ExternalFieldFromData(){
        delete _spl;
    }

    double operator()(const double& t) const override{
        return (*_spl)(t);
    }
private:
    spline_t* _spl;
};

#endif