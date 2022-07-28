#ifndef GRAPHENEMODEL2_HPP
#define GRAPHENEMODEL2_HPP

//This is the model from 
//Faisal, Front. Chem. 10:859405 (2022). doi: 10.3389/fchem.2022.859405

//#include "utils.hpp"
#include "external_field.hpp"

#include <boost/math/differentiation/autodiff.hpp>
using namespace boost::math::differentiation;
using namespace boost::math::differentiation::detail;

class NewField:
public ExternalFieldFromData{
public:
    NewField(const ExternalFieldFromData& ef):
    ExternalFieldFromData(ef){}

    template<typename RealType, size_t Order>
    fvar<RealType,Order> operator()(fvar<RealType, Order> const& cr) const{
        using root_type=typename fvar<RealType,Order>::root_type;
        constexpr size_t order=fvar<RealType,Order>::order_sum;

        //zero order derivative
        root_type const d0=ExternalFieldFromData::operator()(static_cast<root_type>(cr));

        if constexpr (order==0){
            return fvar<RealType,Order>(d0);
        }
        else{
            root_type derivatives[order]{0.0};
            derivatives[0]=d0;
            derivatives[1]=ExternalFieldFromData::derivative(static_cast<root_type>(cr));
            return cr.apply_derivatives(order,[&derivatives](size_t i){return derivatives[i];});
        }
    }

    double operator()(const double& t) const{
        return NewField::operator()(make_fvar<double,0>(t)).derivative(0);
        //return 0.;
    }
};

class GrapheneModel2{
    typedef matrix<complex_t> state_type;
public:
    GrapheneModel2(const double& a,
                   const double& e2p, const double& gamma, const double& s,
                   ExternalField* Ex, ExternalField* Ey,
                   ExternalField* Ax, ExternalField* Ay):
    _a(a),
    _e2p(e2p),_g0(-gamma),_s0(s){
        _Ex=new NewField(*static_cast<ExternalFieldFromData*>(Ex));
        _Ey=new NewField(*static_cast<ExternalFieldFromData*>(Ey));
        _Ax=new NewField(*static_cast<ExternalFieldFromData*>(Ax));
        _Ay=new NewField(*static_cast<ExternalFieldFromData*>(Ay));
    }

    template<typename T>
    T h1(const T& kx, const T& ky) const{
        //return 2.*cos(0.5*kx*_a/sqrt(3.))*cos(0.5*ky*_a);
        return cos(kx*_a/sqrt(3.))+2.*cos(-0.5*kx*_a/sqrt(3.))*cos(0.5*ky*_a);
    }

    template<typename T>
    T h2(const T& kx, const T& ky) const{
        return sin(kx*_a/sqrt(3.))+2.*sin(-0.5*kx*_a/sqrt(3.))*cos(0.5*ky*_a);
    }

    template<typename T>
    T chi(const T& kx, const T& ky) const{
        return atan(h2(kx,ky)/h1(kx,ky));
    }

    template<typename T1, typename T2, typename T3>
    promote<T1,T2,T3> chi_t(const T1& kx0, const T2& ky0, const T3& t) const{
        auto kxt=kx0+(*_Ax)(t);
        auto kyt=ky0+(*_Ay)(t);
        return chi(kxt,kyt);
    }

    double dchidt(const double& kx0, const double& ky0, const double& t) const{
        auto const S=make_fvar<double,1>(t);
        auto const autodiff=chi_t(kx0,ky0,S);
        return autodiff.derivative(1);
    }

    double f(const double& kx, const double& ky) const{
        return sqrt(h1(kx,ky)*h1(kx,ky)+h2(kx,ky)*h2(kx,ky));
    }

    double E1(const double& kx, const double& ky) const{
        return (_e2p-_g0*f(kx,ky))/(1.+_s0*f(kx,ky));
    }

    double E2(const double& kx, const double& ky) const{
        return (_e2p+_g0*f(kx,ky))/(1.-_s0*f(kx,ky));
    }

    void propagate(const state_type &c, state_type &dcdt, const double t,
                   const double& kx0, const double& ky0) const{
        double kxt=kx0+(*_Ax)(t);
        double kyt=ky0+(*_Ay)(t);

        double dchi=dchidt(kx0,ky0,t);

        double V11=(_e2p-_g0*f(kxt,kyt))/(1.+_s0*f(kxt,kyt))-0.5*dchi;
        double V12=0.5*dchi;
        double V21=0.5*dchi;
        double V22=(_e2p+_g0*f(kxt,kyt))/(1.-_s0*f(kxt,kyt))-0.5*dchi;

        dcdt(0,0)=-I*(V11*c(0,0)+V12*c(1,1));
        dcdt(1,1)=-I*(V21*c(0,0)+V22*c(1,1));

        //std::cout<<dchi<<std::endl;

        return;
    }
private:
    double _e2p;
    double _g0;
    double _s0;
    double _a;

    NewField* _Ex;
    NewField* _Ey;
    NewField* _Ax;
    NewField* _Ay;
};

#endif