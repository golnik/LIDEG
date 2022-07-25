#ifndef GRAPHENEMODEL_HPP
#define GRAPHENEMODEL_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include "utils.hpp"
#include "external_field.hpp"

class GrapheneModel{
    friend class WFs;
    typedef matrix<complex_t> state_type;
public:
    GrapheneModel(const double& a,
                  const double& e2p, const double& gamma, const double& s,
                  ExternalField* Ex, ExternalField* Ey):
    _a(a),
    _e2p(e2p),_gamma(gamma),_s(s),
    _Ex{Ex},_Ey{Ey}{}

    ~GrapheneModel(){}

    void get_Dirac_points(std::vector<double>& Dirac_Kx, 
                          std::vector<double>& Dirac_Ky,
                          std::vector<int>&    Dirac_type) const{
        double Kx=2.*M_PI/(sqrt(3.)*_a);
        double Ky=2.*M_PI/(3.*_a);

        //add K point
        Dirac_Kx.push_back(Kx);
        Dirac_Ky.push_back(Ky);
        Dirac_type.push_back(0);

        //add K' point
        Dirac_Kx.push_back(Kx);
        Dirac_Ky.push_back(-Ky);
        Dirac_type.push_back(1);
    }

    complex_t f(const double& kx, const double& ky) const{
        return exp(I*kx*_a/sqrt(3.))+2.*exp(-I*0.5*kx*_a/sqrt(3.))*cos(0.5*ky*_a);
    }

    double phi(const double& kx, const double& ky) const{
        return arg(f(kx,ky));
    }

    double ep(const double& kx, const double& ky) const{
        return (_e2p+_gamma*abs(f(kx,ky)))/(1.+_s*abs(f(kx,ky)));
    }

    double em(const double& kx, const double& ky) const{
        return (_e2p-_gamma*abs(f(kx,ky)))/(1.-_s*abs(f(kx,ky)));
    }

    double d_prefac(const double& kx, const double& ky) const{
        return _a/pow(abs(f(kx,ky)),2.)
               *sqrt((1.-_s*abs(f(kx,ky)))/(1.+_s*abs(f(kx,ky))));
    }

    double dx(const double& kx, const double& ky) const{
        return 1./(2.*sqrt(3.))*d_prefac(kx,ky)
               *(cos(0.5*sqrt(3.)*_a*kx)*cos(0.5*_a*ky)-cos(_a*ky));
    }

    double dy(const double& kx, const double& ky) const{
        return 0.5*d_prefac(kx,ky)
               *sin(0.5*sqrt(3.)*_a*kx)*sin(0.5*_a*ky);
    }

    double px_vv(const double& kx, const double& ky) const{
        return (_a*_gamma/abs(f(kx,ky)))
               *sqrt(3.)*sin(0.5*sqrt(3.)*_a*kx)*cos(0.5*_a*ky);
    }

    double py_vv(const double& kx, const double& ky) const{
        return (_a*_gamma/abs(f(kx,ky)))
                *(cos(0.5*sqrt(3.)*_a*kx)*sin(0.5*_a*ky)+sin(_a*ky));
    }

    double px_cc(const double& kx, const double& ky) const{
        return -px_vv(kx,ky);
    }

    double py_cc(const double& kx, const double& ky) const{
        return -py_vv(kx,ky);
    }

    complex_t px_cv(const double& kx, const double& ky) const{
        return I*(ep(kx,ky)-em(kx,ky))*dx(kx,ky);
    }

    complex_t py_cv(const double& kx, const double& ky) const{
        return I*(ep(kx,ky)-em(kx,ky))*dy(kx,ky);
    }    

    void print(){
        std::cout<<"a: "<<_a<<std::endl;
        std::cout<<"e2p: "<<_e2p<<std::endl;
        std::cout<<"gamma: "<<_gamma<<std::endl;
        std::cout<<"s: "<<_s<<std::endl;
    }

    void propagate(const state_type& rho, state_type& drhodt, const double t,
                    const double& kx_t, const double& ky_t) const{
                        
        double emn[2]={ep(kx_t,ky_t),em(kx_t,ky_t)};

        matrix<complex_t> d_x(2,2);
        d_x(0,0)=0.;
        d_x(1,1)=0.;
        d_x(0,1)=dx(kx_t,ky_t);
        d_x(1,0)=dx(kx_t,ky_t);

        matrix<complex_t> d_y(2,2);
        d_y(0,0)=0.;
        d_y(1,1)=0.;
        d_y(0,1)=dy(kx_t,ky_t);
        d_y(1,0)=dy(kx_t,ky_t);

        //auto comm_x=prod(d_x,rho)-prod(rho,d_x);
        //auto comm_y=prod(d_y,rho)-prod(rho,d_y);

        for(size_t n=0; n<2; n++){
            for(size_t m=0; m<2; m++){
                complex_t res_x=0.;
                complex_t res_y=0.;
                for(size_t mp=0; mp<2; mp++){
                    res_x+=d_x(mp,n)*rho(mp,m)-d_x(m,mp)*rho(n,mp);
                    res_y+=d_y(mp,n)*rho(mp,m)-d_y(m,mp)*rho(n,mp);
                }

                drhodt(n,m)=-I*(
                    (emn[m]-emn[n])*rho(n,m)
                    //-( ((*_Ex)(t)*d_x(m,n)+(*_Ey)(t)*d_y(m,n))*rho(m,m)
                    //  -((*_Ex)(t)*std::conj(d_x(m,n))+(*_Ey)(t)*std::conj(d_y(m,n)))*rho(n,n)
                    // )

                    -(*_Ex)(t)*res_x-(*_Ey)(t)*res_y

                    //-(*_Ex)(t)*comm_x(n,m)//-(*_Ey)(t)*comm_y(m,n)

                );
            }
        }
    }
private:
    double _a;
    double _e2p;
    double _gamma;
    double _s;

    ExternalField* _Ex;
    ExternalField* _Ey;
};

#endif