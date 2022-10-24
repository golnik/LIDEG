#include <iostream>
#include <complex>

#include "external_field.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/math/differentiation/autodiff.hpp>
using namespace boost::math::differentiation;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

using Eigen::MatrixXcd;

class HexagonalTBModel{
    friend class NGraphene;
public:
    HexagonalTBModel(const double& a):
    _a(a){}

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

    complex_t df_dx(const double& kx, const double& ky) const{
        auto const kx_ad=make_fvar<double,1>(kx);
        auto const f_re_ad=f_re(kx_ad,ky);
        auto const f_im_ad=f_im(kx_ad,ky);

        return f_re_ad.derivative(1)+I*f_im_ad.derivative(1);
    }

    complex_t df_dy(const double& kx, const double& ky) const{
        auto const ky_ad=make_fvar<double,1>(ky);
        auto const f_re_ad=f_re(kx,ky_ad);
        auto const f_im_ad=f_im(kx,ky_ad);

        return f_re_ad.derivative(1)+I*f_im_ad.derivative(1);
    }    
private:
    double _a;
};

class NGraphene{
    typedef MatrixXcd matrix_t;
    typedef Eigen::Matrix2cd matrix2D_t;
    typedef Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXcd> solver_t;

    typedef matrix<complex_t> state_type;
public:
    NGraphene(HexagonalTBModel* tbm, const size_t& N,
              double* eps,
              double* g,
              double* s,
              ExternalField* Ex, ExternalField* Ey):
    _tbm{tbm},_N(N),
    _eps{eps},_g{g},_s{s},
    _Ex{Ex},_Ey{Ey}{
        _H0=new matrix_t(2*N,2*N);
        _S=new matrix_t(2*N,2*N);

        _dH0_dx=new matrix_t(2*N,2*N);
        _dH0_dy=new matrix_t(2*N,2*N);

        _solver=new solver_t;
    }

    int solve(const double& kx, const double& ky) const{
        double epsA=_eps[0];
        double epsB=_eps[1];

        double g0=_g[0];
        double g1=_g[1];

        double s0=_s[0];
        double s1=_s[1];

        matrix2D_t Hii;
        matrix2D_t Hij;

        matrix2D_t dHii_dx;
        matrix2D_t dHij_dx;
        matrix2D_t dHii_dy;
        matrix2D_t dHij_dy;

        matrix2D_t Sii;
        matrix2D_t Sij;

        complex_t fval =_tbm->f(kx,ky);
        complex_t df_dx=_tbm->df_dx(kx,ky);
        complex_t df_dy=_tbm->df_dy(kx,ky);

        Hii(0,0)=epsA;
        Hii(0,1)=g0*fval;
        Hii(1,0)=g0*std::conj(fval);
        Hii(1,1)=epsB;

        Hij(1,0)=g1;

        dHii_dx(0,1)=g0*df_dx;
        dHii_dx(1,0)=g0*std::conj(df_dx);

        dHii_dy(0,1)=g0*df_dy;
        dHii_dy(1,0)=g0*std::conj(df_dy);

        Sii(0,0)=1.;
        Sii(0,1)=s0*fval;
        Sii(1,0)=s0*std::conj(fval);
        Sii(1,1)=1.;

        Sij(1,0)=s1;

        for(size_t i=0; i<_N; i++){
            _H0->block<2,2>(2*i,2*i)=Hii;
            _S->block<2,2>(2*i,2*i)=Sii;

            _dH0_dx->block<2,2>(2*i,2*i)=dHii_dx;
            _dH0_dy->block<2,2>(2*i,2*i)=dHii_dy;

            if(i<(_N-1)){
                _H0->block<2,2>(2*i,2*(i+1))=Hij;
                _H0->block<2,2>(2*(i+1),2*i)=Hij.adjoint();

                _S->block<2,2>(2*i,2*(i+1))=Sij;
                _S->block<2,2>(2*(i+1),2*i)=Sij.adjoint();

                //_dH0->block<2,2>(2*i,2*(i+1))=dHij;
                //_dH0->block<2,2>(2*(i+1),2*i)=dHij.adjoint();
            }
        }

        //std::cout<<(*_H0)<<std::endl;
        //std::cout<<(*_S)<<std::endl;
        //std::cout<<(*_dH0)<<std::endl;

        _solver->compute(*_H0,*_S);

        return 0;
    }
    
    std::vector<double> getEnergies() const{
        auto evals=_solver->eigenvalues();
        std::vector<double> res(evals.data(),evals.data()+evals.rows()*evals.cols());
        return res;
    }

    double getDipole(const size_t& ist, const size_t& jst) const{
        if(ist==jst){
            return 0.;
        }
        else{
            auto evals=_solver->eigenvalues();
            auto evecs=_solver->eigenvectors();

            double ei=evals(ist);
            double ej=evals(jst);

            double dE=fabs(ei-ej);

            if(dE>=1.e-7){
                auto veci=evecs.col(ist);
                auto vecj=evecs.col(jst);

                //complex_t dip=I*veci.dot(_S->inverse()*(*_dH0_dx)*vecj)/(ei-ej);
                complex_t dip=I*veci.dot((*_dH0_dx)*vecj)/(ei-ej);

                return std::real(dip);
            }
            else{
                return 1.;
            }
        }
    }

    //std::vector<double> 

    void propagate(const state_type& rho, state_type& drhodt, const double t,
                const double& kx_t, const double& ky_t) const{
        size_t Nst=2*_N;

        //this->solve(kx_t,ky_t);
        auto emn=this->getEnergies();

        matrix<complex_t> d_x(Nst,Nst);
        matrix<complex_t> d_y(Nst,Nst);

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=0; jst<Nst; jst++){
                d_x(ist,jst)=this->getDipole(ist,jst);
                d_y(ist,jst)=d_x(ist,jst);
            }
        }

        //auto comm_x=prod(d_x,rho)-prod(rho,d_x);
        //auto comm_y=prod(d_y,rho)-prod(rho,d_y);

        for(size_t n=0; n<Nst; n++){
            for(size_t m=0; m<Nst; m++){
                complex_t res_x=0.;
                complex_t res_y=0.;
                for(size_t mp=0; mp<Nst; mp++){
                    res_x+=d_x(mp,n)*rho(mp,m)-d_x(m,mp)*rho(n,mp);
                    res_y+=d_y(mp,n)*rho(mp,m)-d_y(m,mp)*rho(n,mp);
                }

                //std::cout<<res_x<<" "<<res_y<<std::endl;

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
        return;
    }
private:
    HexagonalTBModel* _tbm;
    size_t _N;

    double* _eps;
    double* _g;
    double* _s;

    matrix_t* _H0;
    matrix_t* _S;
    
    matrix_t* _dH0_dx;
    matrix_t* _dH0_dy;

    solver_t* _solver;

    ExternalField* _Ex;
    ExternalField* _Ey;    
};