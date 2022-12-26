#ifndef GRAPHENEMODEL_HPP
#define GRAPHENEMODEL_HPP

#include <iostream>
#include <cmath>
#include <complex>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/math/differentiation/autodiff.hpp>
using namespace boost::math::differentiation;

//#include "utils.hpp"
#include "utils/grid.hpp"
#include "external_field.hpp"

#include "graphenemodel.hpp"

class GrapheneModel:
public Graphene{
    friend class WFs;
    friend class WFs_grid;
    typedef matrix<complex_t> state_type;
public:
    GrapheneModel(const double& a,
                  const double& e2p, const double& gamma, const double& s,
                  const double& Td,
                  ExternalField* Ex, ExternalField* Ey):
    _a(a),
    _e2p(e2p),_gamma(gamma),_s(s),
    _Td(Td),
    _Ex{Ex},_Ey{Ey}{}

    virtual ~GrapheneModel(){}

    /*complex_t f(const double& kx, const double& ky) const{
        return exp(I*kx*_a/sqrt(3.))+2.*exp(-I*0.5*kx*_a/sqrt(3.))*cos(0.5*ky*_a);
    }*/

    template<typename T1, typename T2>
    promote<T1,T2> f_re(const T1& kx, const T2& ky) const{
        return cos(kx*_a/sqrt(3.))+2.*cos(-0.5*kx*_a/sqrt(3.))*cos(0.5*ky*_a);
    }

    template<typename T1, typename T2>
    promote<T1,T2> f_im(const T1& kx, const T2& ky) const{
        return sin(kx*_a/sqrt(3.))+2.*sin(-0.5*kx*_a/sqrt(3.))*cos(0.5*ky*_a);
    }

    complex_t f(const double& kx, const double& ky) const{
        return f_re(kx,ky)+I*f_im(kx,ky);
    }

    template<typename T1, typename T2>
    promote<T1,T2> fabs(const T1& kx, const T2& ky) const{
        return sqrt(f_re(kx,ky)*f_re(kx,ky)+f_im(kx,ky)*f_im(kx,ky));
    }    

    double phi(const double& kx, const double& ky) const{
        return atan2(f_im(kx,ky),f_re(kx,ky)); //arg(f(kx,ky));
    }

    template<typename T1, typename T2>
    promote<T1,T2> ep(const T1& kx, const T2& ky) const{
        return (_e2p+_gamma*fabs(kx,ky))/(1.+_s*fabs(kx,ky));
    }

    template<typename T1, typename T2>
    promote<T1,T2> em(const T1& kx, const T2& ky) const{
        return (_e2p-_gamma*fabs(kx,ky))/(1.-_s*fabs(kx,ky));
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
        //return (_a*_gamma/abs(f(kx,ky)))
        //       *sqrt(3.)*sin(0.5*sqrt(3.)*_a*kx)*cos(0.5*_a*ky);
        auto const kx_ad=make_fvar<double,1>(kx);
        auto const autodiff=ep(kx_ad,ky);
        return autodiff.derivative(1);
    }

    double py_vv(const double& kx, const double& ky) const{
//        return (_a*_gamma/abs(f(kx,ky)))
//                *(cos(0.5*sqrt(3.)*_a*kx)*sin(0.5*_a*ky)+sin(_a*ky));
        auto const ky_ad=make_fvar<double,1>(ky);
        auto const autodiff=ep(kx,ky_ad);
        return autodiff.derivative(1);
    }

    double px_cc(const double& kx, const double& ky) const{
        //return -px_vv(kx,ky);
        auto const kx_ad=make_fvar<double,1>(kx);
        auto const autodiff=em(kx_ad,ky);
        return autodiff.derivative(1);
    }

    double py_cc(const double& kx, const double& ky) const{
        //return -py_vv(kx,ky);
        auto const ky_ad=make_fvar<double,1>(ky);
        auto const autodiff=em(kx,ky_ad);
        return autodiff.derivative(1);
    }

    double get_energy(const double& kx, const double& ky, const size_t& ist) const override{
        double res=0.;
        if(ist==0){
            res=ep(kx,ky);
        }
        else if(ist==1){
            res=em(kx,ky);
        }
        else{
            throw std::string("Only two energy levels are available in this model!");
        }
        return res;
    }

    double get_dipole(const double& kx, const double& ky, const size_t& ist, const size_t& jst, const size_t& dir) const override{
        double res=0.;
        if(ist!=jst){
            if(dir==0){
                res=dx(kx,ky);
            }
            else if(dir==1){
                res=dy(kx,ky);
            }
        }
        return res;
    }

    double get_energy_grad(const double& kx, const double& ky, const size_t& ist, const size_t& dir) const override{
        double res=0.;
        if(dir==0){
            if(ist==0){
                res=px_vv(kx,ky);
            }
            else if(ist==1){
                res=px_cc(kx,ky);
            }
        }
        else if(dir==1){
            if(ist==0){
                res=py_vv(kx,ky);
            }
            else if(ist==1){
                res=py_cc(kx,ky);
            }
        }
        return res;
    }

    complex_t px_cv(const double& kx, const double& ky) const{
        return I*(em(kx,ky)-ep(kx,ky))*dx(kx,ky);
    }

    complex_t py_cv(const double& kx, const double& ky) const{
        return I*(em(kx,ky)-ep(kx,ky))*dy(kx,ky);
    }    

    void print(){
        std::cout<<"a: "<<_a<<std::endl;
        std::cout<<"e2p: "<<_e2p<<std::endl;
        std::cout<<"gamma: "<<_gamma<<std::endl;
        std::cout<<"s: "<<_s<<std::endl;
    }

    void write_energies_to_file(const std::string& fname, Grid2D* kxygrid) const{
        std::ofstream out(fname);

        size_t Nkx=kxygrid->size1();
        size_t Nky=kxygrid->size2();

        out<<std::scientific;
        out<<std::setprecision(8);
        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                double kx=(*kxygrid)(ikx,iky)[0];
                double ky=(*kxygrid)(ikx,iky)[1];

                out<<std::setw(20)<<kx<<std::setw(20)<<ky;
                out<<std::setw(20)<<ep(kx,ky);
                out<<std::setw(20)<<em(kx,ky);
                out<<std::endl;
            }
        }

        out.close();
        return;
    }

    void write_dipoles_to_file(const std::string& fname, Grid2D* kxygrid) const{
        std::ofstream out(fname);

        size_t Nkx=kxygrid->size1();
        size_t Nky=kxygrid->size2();

        out<<std::scientific;
        out<<std::setprecision(8);
        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                double kx=(*kxygrid)(ikx,iky)[0];
                double ky=(*kxygrid)(ikx,iky)[1];

                out<<std::setw(20)<<kx<<std::setw(20)<<ky;
                out<<std::setw(20)<<std::abs(dx(kx,ky));
                out<<std::setw(20)<<std::abs(dy(kx,ky));
                out<<std::endl;
            }
        }

        out.close();
        return;
    }    

    void propagate(const state_type& rho, state_type& drhodt, const double t,
                    const double& kx_t, const double& ky_t) const override{
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

                if(m!=n) drhodt(n,m)+=-(1./_Td)*rho(n,m);
            }
        }
    }

    size_t nstates() const override{
        return 2;
    }
private:
    double _a;
    double _e2p;
    double _gamma;
    double _s;
    double _Td;

    ExternalField* _Ex;
    ExternalField* _Ey;
};

#endif