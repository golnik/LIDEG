#include <iostream>
#include <complex>
#include <memory>

#include "external_field.hpp"

#include "utils/grid.hpp"

#include "graphenemodel.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/math/differentiation/autodiff.hpp>
using namespace boost::math::differentiation;

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

using Eigen::MatrixXcd;

#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>

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

        //add K point
        /*Dirac_Kx.push_back(-Kx);
        Dirac_Ky.push_back(Ky);
        Dirac_type.push_back(0);

        //add K' point
        Dirac_Kx.push_back(-Kx);
        Dirac_Ky.push_back(-Ky);
        Dirac_type.push_back(1);*/

        //add K point
        /*Dirac_Kx.push_back(4.*M_PI/(3.*_a));
        Dirac_Ky.push_back(0.);
        Dirac_type.push_back(0);

        //add' K point
        Dirac_Kx.push_back(-4.*M_PI/(3.*_a));
        Dirac_Ky.push_back(0.);
        Dirac_type.push_back(0);*/
    }    
private:
    double _a;
};

class NGraphene:
public Graphene{
    typedef MatrixXcd matrix_t;
    typedef Eigen::Matrix2cd matrix2D_t;
    typedef Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXcd> solver_t;

    typedef matrix<complex_t> state_type;
public:
    NGraphene(HexagonalTBModel* tbm, const size_t& N,
              Grid2D* kxygrid,
              double* eps,
              double* g,
              double* s,
              ExternalField* Ex, ExternalField* Ey):
    _tbm{tbm},_N(N),
    _eps{eps},_g{g},_s{s},
    _Ex{Ex},_Ey{Ey}{
        size_t Nst=2*N;

        matrix_t zero_matrix(Nst,Nst);
        zero_matrix=matrix_t::Zero(Nst,Nst);

        _H0=new matrix_t(2*N,2*N);
        *_H0=zero_matrix;
        _S=new matrix_t(2*N,2*N);
        *_S=zero_matrix;

        _dH0.push_back(new matrix_t(2*N,2*N));
        _dH0.push_back(new matrix_t(2*N,2*N));
        *_dH0[0]=zero_matrix;
        *_dH0[1]=zero_matrix;

        _solver=new solver_t;

        //fit energies and dipoles by splines
        size_t Nkx=kxygrid->size1();
        size_t Nky=kxygrid->size2();

        double xmin=(*kxygrid)(0,0)[0];
        double ymin=(*kxygrid)(0,Nky-1)[1];
        double xmax=(*kxygrid)(Nkx-1,Nky-1)[0];
        double ymax=(*kxygrid)(Nkx-1,0)[1];

        auto kx_grid=create_grid(xmin,xmax,Nkx,0);
        auto ky_grid=create_grid(ymin,ymax,Nky,0);

        SPLINTER::DenseVector kxky(2);
        std::vector<SPLINTER::DataTable> e_data(Nst);
        int Ndip=int(Nst*(Nst-1)/2);
        std::vector<SPLINTER::DataTable> d_data_x(Ndip);
        std::vector<SPLINTER::DataTable> d_data_y(Ndip);

        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                double kx=kx_grid[ikx];
                double ky=ky_grid[iky];

                this->solve(kx,ky);
                auto evals=_solver->eigenvalues();
                auto evecs=_solver->eigenvalues();

                kxky(0)=kx;
                kxky(1)=ky;

                //std::cout<<kx<<" "<<ky<<" ";

                //std::cout<<"evals: "<<std::endl;
                //std::cout<<evals<<std::endl;
                //std::cout<<"evecs: "<<std::endl;
                //std::cout<<evecs<<std::endl;

                size_t indx=0;
                for(size_t ist=0; ist<Nst; ist++){
                    e_data[ist].addSample(kxky,evals[ist]);
                    for(size_t jst=ist+1; jst<Nst; jst++){
                        double dip_x=this->computeDipole(ist,jst,0);
                        double dip_y=this->computeDipole(ist,jst,1);

                        d_data_x[indx].addSample(kxky,dip_x);
                        d_data_y[indx].addSample(kxky,dip_y);

                        //std::cout<<dip_x<<" ";
                        //std::cout<<dip_y<<" ";

                        indx++;
                    }
                }
                //std::cout<<std::endl;
            }
        }

        //fit the states
        size_t indx=0;
        for(size_t ist=0; ist<Nst; ist++){
            std::cout<<"Interpolation of state "<<ist+1<<" by BSplines"<<std::endl;
            e_spl.push_back(std::make_shared<SPLINTER::BSpline>(SPLINTER::BSpline::Builder(e_data[ist]).degree(1).build()));
            for(size_t jst=ist+1; jst<Nst; jst++){
                d_spl_x.push_back(std::make_shared<SPLINTER::BSpline>(SPLINTER::BSpline::Builder(d_data_x[indx]).degree(1).build()));
                d_spl_y.push_back(std::make_shared<SPLINTER::BSpline>(SPLINTER::BSpline::Builder(d_data_y[indx]).degree(1).build()));
                indx++;
            }            
        }

        /*double kx=0.;
        double ky=0.0;
        kxky(0)=kx;
        kxky(1)=ky;
        for(size_t ist=0; ist<Nstates; ist++){
            //double res=e_spl[ist].eval(kxky);//this->spl_evaluate(e_spl[ist],kx,ky);
            double res=this->spl_evaluate(*e_spl[ist],kx,ky);
            std::cout<<ist<<" "<<res<<std::endl;
        }*/
    }

    void propagate(const state_type& rho, state_type& drhodt, const double t,
                const double& kx_t, const double& ky_t) const override{
        size_t Nst=2*_N;

        std::vector<double> emn(Nst);
        for(size_t ist=0; ist<Nst; ist++){
            emn[ist]=this->spl_evaluate(*e_spl[ist],kx_t,ky_t);
        }

        //this->solve(kx_t,ky_t);
        //auto emn=this->getEnergies();

        matrix<complex_t> d_x(Nst,Nst);
        matrix<complex_t> d_y(Nst,Nst);

        size_t indx=0;
        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=ist+1; jst<Nst; jst++){
                d_x(ist,jst)=this->spl_evaluate(*d_spl_x[indx],kx_t,ky_t);
                d_x(jst,ist)=d_x(ist,jst);

                d_y(ist,jst)=this->spl_evaluate(*d_spl_y[indx],kx_t,ky_t);
                d_y(jst,ist)=d_y(ist,jst);

                indx++;
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

    size_t nstates() const override{
        return 2*_N;
    }

    void write_energies_to_file(const std::string& fname, Grid2D* kxygrid) const{
        std::ofstream out(fname);

        size_t Nkx=kxygrid->size1();
        size_t Nky=kxygrid->size2();
        size_t Nst=2*_N;

        out<<std::scientific;
        out<<std::setprecision(8);
        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                double kx=(*kxygrid)(ikx,iky)[0];
                double ky=(*kxygrid)(ikx,iky)[1];

                out<<std::setw(20)<<kx<<std::setw(20)<<ky;
                for(size_t ist=0; ist<Nst; ist++){
                    double en=this->spl_evaluate(*e_spl[ist],kx,ky);
                    out<<std::setw(20)<<en;
                }
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
        size_t Nst=2*_N;

        out<<std::scientific;
        out<<std::setprecision(8);
        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                double kx=(*kxygrid)(ikx,iky)[0];
                double ky=(*kxygrid)(ikx,iky)[1];

                this->solve(kx,ky);

                out<<std::setw(20)<<kx<<std::setw(20)<<ky;
                size_t indx=0;
                for(size_t ist=0; ist<Nst; ist++){
                    for(size_t jst=ist+1; jst<Nst; jst++){
                        //double dip_x=this->spl_evaluate(*d_spl_x[indx],kx,ky);
                        //double dip_y=this->spl_evaluate(*d_spl_y[indx],kx,ky);
                        double dip_x=this->computeDipole(ist,jst,0);
                        double dip_y=this->computeDipole(ist,jst,1);
                        out<<std::setw(20)<<dip_x;
                        out<<std::setw(20)<<dip_y;
                        indx++;
                    }
                }
                out<<std::endl;
            }
        }

        out.close();
        return;
    }
//private:
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

            _dH0[0]->block<2,2>(2*i,2*i)=dHii_dx;
            _dH0[1]->block<2,2>(2*i,2*i)=dHii_dy;

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
    
    /*std::vector<double> getEnergies() const{
        auto evals=_solver->eigenvalues();
        std::vector<double> res(evals.data(),evals.data()+evals.rows()*evals.cols());
        return res;
    }*/

    double computeDipole(const size_t& ist, const size_t& jst, const size_t& dir) const{
        double res=0.;
        if(ist!=jst){
            auto evals=_solver->eigenvalues();
            auto evecs=_solver->eigenvectors();

            double ei=evals(ist);
            double ej=evals(jst);

            auto veci=evecs.col(ist);
            auto vecj=evecs.col(jst);

            //complex_t dip=I*veci.dot(_S->inverse()*(*_dH0[dir])*vecj)/(ei-ej);
            //complex_t dip=I*veci.dot((*_dH0[dir])*vecj)/(ei-ej);
            //res=std::real(dip);
            
            //std::cout<<"S_inverse: "<<std::endl;
            //std::cout<<_S->inverse()<<std::endl;

            //std::cout<<"dH0: "<<std::endl;
            //std::cout<<*_dH0[dir]<<std::endl;

            complex_t tmp=I*veci.dot(_S->inverse()*(*_dH0[dir])*vecj);
            res=-std::abs(tmp)/(ei-ej);
        }
        return res;
    }

    //std::vector<double> 

    double spl_evaluate(const SPLINTER::BSpline& spl, const double& kx, const double& ky) const{
        SPLINTER::DenseVector kxky(2);
        kxky(0)=kx;
        kxky(1)=ky;
        return spl.eval(kxky);
    }

    HexagonalTBModel* _tbm;
    size_t _N;

    double* _eps;
    double* _g;
    double* _s;

    matrix_t* _H0;
    matrix_t* _S;
    
    std::vector<matrix_t*> _dH0;

    solver_t* _solver;

    ExternalField* _Ex;
    ExternalField* _Ey;

    std::vector<std::shared_ptr<SPLINTER::BSpline>> e_spl;
    std::vector<std::shared_ptr<SPLINTER::BSpline>> d_spl_x;
    std::vector<std::shared_ptr<SPLINTER::BSpline>> d_spl_y;
};