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
    typedef Eigen::MatrixXcd matrix_t;
    typedef Eigen::Matrix2cd matrix2D_t;
    typedef Eigen::VectorXd Evector_t;

    typedef matrix<complex_t> state_type;
public:
    NGraphene(HexagonalTBModel* tbm, const size_t& N,
              Grid2D* kxygrid,
              const std::vector<double>& eps,
              const std::vector<double>& g,
              const std::vector<double>& s,
              const double& Td,
              ExternalField* Ex, ExternalField* Ey):
    _tbm{tbm},_N(N),
    _eps(eps),_g(g),_s(s),
    _Td(Td),
    _Ex{Ex},_Ey{Ey}{
        /*size_t Nst=2*N;

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

                Evector_t evals;
                matrix_t evecs;
                this->solve(evals,evecs,kx,ky);

                matrix_t Dx=matrix_t::Zero(Nst,Nst);
                matrix_t Dy=matrix_t::Zero(Nst,Nst);
                this->computeDipole(Dx,Dy,kx,ky);

                kxky(0)=kx;
                kxky(1)=ky;

                size_t indx=0;
                for(size_t ist=0; ist<Nst; ist++){
                    double en=evals[ist];
                    e_data[ist].addSample(kxky,en);
                    for(size_t jst=ist+1; jst<Nst; jst++){
                        d_data_x[indx].addSample(kxky,std::real(Dx(ist,jst)));
                        d_data_y[indx].addSample(kxky,std::real(Dy(ist,jst)));
                        indx++;
                    }
                }
            }
        }

        mapping=new Eigen::MatrixXi(Nst,Nst);

        //fit the states
        size_t indx=0;
        for(size_t ist=0; ist<Nst; ist++){
            std::cout<<"Interpolation of state "<<ist+1<<" by BSplines"<<std::endl;
            e_spl.push_back(std::make_shared<SPLINTER::BSpline>(SPLINTER::BSpline::Builder(e_data[ist]).degree(1).build()));
            for(size_t jst=ist+1; jst<Nst; jst++){
                d_spl_x.push_back(std::make_shared<SPLINTER::BSpline>(SPLINTER::BSpline::Builder(d_data_x[indx]).degree(1).build()));
                d_spl_y.push_back(std::make_shared<SPLINTER::BSpline>(SPLINTER::BSpline::Builder(d_data_y[indx]).degree(1).build()));

                (*mapping)(ist,jst)=indx;

                indx++;
            }            
        }*/
    }

    double get_energy(const double& kx, const double& ky, const size_t& ist) const override{
        Evector_t evals;
        matrix_t evecs;
        this->solve(evals,evecs,kx,ky);
        return evals[ist];
        //return this->spl_evaluate(*e_spl[ist],kx,ky);
    }

    double get_dipole(const double& kx, const double& ky, const size_t& ist, const size_t& jst, const size_t& dir) const override{
        size_t Nst=2*_N;

        matrix_t Dx=matrix_t::Zero(Nst,Nst);
        matrix_t Dy=matrix_t::Zero(Nst,Nst);
        this->computeDipole(Dx,Dy,kx,ky);

        double res=0.;
        if(dir==0){
            res=std::real(Dx(ist,jst));
        }
        else if(dir==1){
            res=std::real(Dy(ist,jst));
        }

        return res;

        /*size_t indx=(*mapping)(ist,jst);

        double res=0.;
        if(ist!=jst){
            if(dir==0){
                res=this->spl_evaluate(*d_spl_x[indx],kx,ky);
            }
            else if(dir==1){
                res=this->spl_evaluate(*d_spl_y[indx],kx,ky);
            }
        }
        return res;*/
    }

    double get_energy_grad(const double& kx, const double& ky, const size_t& ist, const size_t& dir) const override{
        Evector_t evals;
        matrix_t evecs;
        this->solve(evals,evecs,kx,ky);

        size_t Nst=this->nstates();

        matrix_t H=matrix_t::Zero(Nst,Nst);
        matrix_t S=matrix_t::Zero(Nst,Nst);

        this->compute_matrices(H,S,kx,ky);

        matrix_t dH=matrix_t::Zero(Nst,Nst);
        matrix_t dS=matrix_t::Zero(Nst,Nst);

        this->compute_derivatives(dH,dS,kx,ky,dir);

        auto veci=evecs.col(ist);

        complex_t res=veci.dot(S.inverse()*dH*veci);

        return std::real(res);
    }

    vector_t get_state(const double& kx, const double& ky, const size_t& ist) const override{
        Evector_t evals;
        matrix_t evecs;
        this->solve(evals,evecs,kx,ky);

        size_t Nst=this->nstates();

        vector_t res(Nst);

        auto veci=evecs.col(ist);
        for(size_t j=0; j<this->nstates(); j++){
            res(j)=veci(j);
        }

        return res;
    }

    void propagate(const state_type& rho, state_type& drhodt, const double t,
                const double& kx_t, const double& ky_t) const override{
        //size_t Nst=2*_N;
        size_t Nst=this->nstates();

        Evector_t evals;
        matrix_t evecs;
        matrix_t Dx=matrix_t::Zero(Nst,Nst);
        matrix_t Dy=matrix_t::Zero(Nst,Nst);

        this->solve(evals,evecs,kx_t,ky_t);

        matrix_t H=matrix_t::Zero(Nst,Nst);
        matrix_t S=matrix_t::Zero(Nst,Nst);

        this->compute_matrices(H,S,kx_t,ky_t);

        matrix_t dH_dx=matrix_t::Zero(Nst,Nst);
        matrix_t dS_dx=matrix_t::Zero(Nst,Nst);
        matrix_t dH_dy=matrix_t::Zero(Nst,Nst);
        matrix_t dS_dy=matrix_t::Zero(Nst,Nst);

        this->compute_derivatives(dH_dx,dS_dx,kx_t,ky_t,0);
        this->compute_derivatives(dH_dy,dS_dy,kx_t,ky_t,1);

        matrix_t X=S.inverse()*dH_dx;
        matrix_t Y=S.inverse()*dH_dy;

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=ist+1; jst<Nst; jst++){
                double ei=evals(ist);
                double ej=evals(jst);

                complex_t dip_x=0.;
                complex_t dip_y=0.;

                if(abs(ei-ej)>1.e-10){
                    auto veci=evecs.col(ist);
                    auto vecj=evecs.col(jst);

                    if(std::real(veci(0))<0){
                        veci*=-1;
                    }

                    if(std::real(vecj(0))<0){
                        vecj*=-1;
                    }

                    dip_x=I*veci.dot(X*vecj)/(ei-ej);
                    dip_y=I*veci.dot(Y*vecj)/(ei-ej);
                }

                Dx(ist,jst)=dip_x;
                Dx(jst,ist)=std::conj(dip_x);

                Dy(ist,jst)=dip_y;
                Dy(jst,ist)=std::conj(dip_y);
            }
        }

        std::vector<double> E(Nst);
        for(size_t ist=0; ist<Nst; ist++){
            E[ist]=evals(ist);//this->get_energy(kx_t,ky_t,ist);
        }

        matrix<complex_t> d_x(Nst,Nst);
        matrix<complex_t> d_y(Nst,Nst);

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=0; jst<Nst; jst++){
                d_x(ist,jst)=0.;
                d_y(ist,jst)=0.;
            }
        }

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=ist+1; jst<Nst; jst++){
                double dx=std::real(Dx(ist,jst));//this->get_dipole(kx_t,ky_t,ist,jst,0);
                double dy=std::real(Dy(ist,jst));//this->get_dipole(kx_t,ky_t,ist,jst,1);

                d_x(ist,jst)=dx;
                d_x(jst,ist)=dx;

                d_y(ist,jst)=dy;
                d_y(jst,ist)=dy;
            }
        }

        auto comm_x=prod(d_x,rho)-prod(rho,d_x);
        auto comm_y=prod(d_y,rho)-prod(rho,d_y);

        for(size_t m=0; m<Nst; m++){
            for(size_t n=0; n<Nst; n++){
                drhodt(m,n)=-I*(
                    (E[m]-E[n])*rho(m,n)
                    +(*_Ex)(t)*comm_x(m,n)+(*_Ey)(t)*comm_y(m,n)
                );

                if(m!=n) drhodt(m,n)+=-(1./_Td)*rho(m,n);
            }
        }
        
        return;
    }

    size_t nstates() const override{
        return 2*_N;
    }

private:
    void compute_matrices(matrix_t& H, matrix_t& S, const double& kx, const double& ky) const{
        double epsA=_eps[0];
        double epsB=_eps[1];

        double g0=_g[0];
        double g1=_g[1];
        double g3=_g[2];
        double g4=_g[3];

        double s0=_s[0];
        double s1=_s[1];

        matrix2D_t Hii=matrix2D_t::Zero();
        matrix2D_t Hij=matrix2D_t::Zero();

        matrix2D_t Sii=matrix2D_t::Zero();
        matrix2D_t Sij=matrix2D_t::Zero();

        complex_t fval =_tbm->f(kx,ky);

        Hii(0,0)=epsA;
        Hii(0,1)=g0*fval;
        Hii(1,0)=g0*std::conj(fval);
        Hii(1,1)=epsB;

        Hij(1,0)=g1;
        Hij(0,1)=g3*std::conj(fval);
        Hij(0,0)=g4*fval;
        Hij(1,1)=g4*fval;

        Sii(0,0)=1.;
        Sii(0,1)=s0*fval;
        Sii(1,0)=s0*std::conj(fval);
        Sii(1,1)=1.;

        Sij(1,0)=s1;

        for(size_t i=0; i<_N; i++){
            H.block<2,2>(2*i,2*i)=Hii;
            S.block<2,2>(2*i,2*i)=Sii;

            if(i<(_N-1)){
                H.block<2,2>(2*i,2*(i+1))=Hij;
                H.block<2,2>(2*(i+1),2*i)=Hij.adjoint();

                S.block<2,2>(2*i,2*(i+1))=Sij;
                S.block<2,2>(2*(i+1),2*i)=Sij.adjoint();
            }
        }

        return;
    }

    void compute_derivatives(matrix_t& dH, matrix_t& dS, const double& kx, const double& ky, const size_t& dir) const{
        double epsA=_eps[0];
        double epsB=_eps[1];

        double g0=_g[0];
        double g1=_g[1];
        double g3=_g[2];
        double g4=_g[3];

        double s0=_s[0];
        double s1=_s[1];

        matrix2D_t dHii=matrix2D_t::Zero();
        matrix2D_t dHij=matrix2D_t::Zero();

        complex_t df;
        if(dir==0){
            df=_tbm->df_dx(kx,ky);
        }
        else if(dir==1){
            df=_tbm->df_dy(kx,ky);
        }

        dHii(0,1)=g0*df;
        dHii(1,0)=g0*std::conj(df);

        dHij(0,1)=g3*std::conj(df);
        dHij(0,0)=g4*df;
        dHij(1,1)=g4*df;

        for(size_t i=0; i<_N; i++){
            dH.block<2,2>(2*i,2*i)=dHii;

            if(i<(_N-1)){
                dH.block<2,2>(2*i,2*(i+1))=dHij;
                dH.block<2,2>(2*(i+1),2*i)=dHij.adjoint();
            }
        }

        return;
    }

    int solve(Evector_t& evals, matrix_t& evecs, const double& kx, const double& ky) const{
        size_t Nst=this->nstates();

        matrix_t H=matrix_t::Zero(Nst,Nst);
        matrix_t S=matrix_t::Zero(Nst,Nst);

        this->compute_matrices(H,S,kx,ky);

        Eigen::GeneralizedSelfAdjointEigenSolver<matrix_t> solver;
        solver.compute(H,S);
        
        evals=solver.eigenvalues();
        evecs=solver.eigenvectors();

        return 0;
    }

    void computeDipole(matrix_t& Dx, matrix_t& Dy, const double& kx, const double& ky) const{
        Evector_t evals;
        matrix_t evecs;

        this->solve(evals,evecs,kx,ky);

        size_t Nst=this->nstates();

        matrix_t H=matrix_t::Zero(Nst,Nst);
        matrix_t S=matrix_t::Zero(Nst,Nst);

        this->compute_matrices(H,S,kx,ky);

        matrix_t dH_dx=matrix_t::Zero(Nst,Nst);
        matrix_t dS_dx=matrix_t::Zero(Nst,Nst);
        matrix_t dH_dy=matrix_t::Zero(Nst,Nst);
        matrix_t dS_dy=matrix_t::Zero(Nst,Nst);

        this->compute_derivatives(dH_dx,dS_dx,kx,ky,0);
        this->compute_derivatives(dH_dy,dS_dy,kx,ky,1);

        matrix_t X=S.inverse()*dH_dx;
        matrix_t Y=S.inverse()*dH_dy;

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=ist+1; jst<Nst; jst++){
                double ei=evals(ist);
                double ej=evals(jst);

                complex_t dip_x=0.;
                complex_t dip_y=0.;

                if(abs(ei-ej)>1.e-10){
                    auto veci=evecs.col(ist);
                    auto vecj=evecs.col(jst);

                    if(std::real(veci(0))<0){
                        veci*=-1;
                    }

                    if(std::real(vecj(0))<0){
                        vecj*=-1;
                    }

                    dip_x=I*veci.dot(X*vecj)/(ei-ej);
                    dip_y=I*veci.dot(Y*vecj)/(ei-ej);
                }

                Dx(ist,jst)=dip_x;
                Dx(jst,ist)=std::conj(dip_x);

                Dy(ist,jst)=dip_y;
                Dy(jst,ist)=std::conj(dip_y);
            }
        }

        return;
    }

    /*double spl_evaluate(const SPLINTER::BSpline& spl, const double& kx, const double& ky) const{
        SPLINTER::DenseVector kxky(2);
        kxky(0)=kx;
        kxky(1)=ky;
        return spl.eval(kxky);
    }*/

    HexagonalTBModel* _tbm;
    size_t _N;

    std::vector<double> _eps;
    std::vector<double> _g;
    std::vector<double> _s;
    double _Td;

    ExternalField* _Ex;
    ExternalField* _Ey;

    /*std::vector<std::shared_ptr<SPLINTER::BSpline>> e_spl;
    std::vector<std::shared_ptr<SPLINTER::BSpline>> d_spl_x;
    std::vector<std::shared_ptr<SPLINTER::BSpline>> d_spl_y;

    Eigen::MatrixXi* mapping;*/
};