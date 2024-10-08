#ifndef WFSnew_HPP
#define WFSnew_HPP

#include <cmath>
#include <vector>

#include "utils/utils.hpp"

#include <boost/math/differentiation/autodiff.hpp>
using namespace boost::math::differentiation;

#define THRS 1.e-9

class Orbital{
public:
    Orbital(){}
    ~Orbital(){}

    virtual double p(const double& rx, const double& ry, const double& rz) const=0;

    virtual double dp(const double& rx, const double& ry, const double& rz, const size_t& dir) const{
        return 0;
    }
};

class Pzorb_normal:
public Orbital{
public:
    Pzorb_normal(const double& Z):
    _Z(Z){
        double a=1./_Z;
        _norm=sqrt(1./pow(a,3.))/(4.*a*sqrt(2.*M_PI));
    }

    double p(const double& rx, const double& ry, const double& rz) const override{
        return psi(rx,ry,rz);
    }

    double dp(const double& rx, const double& ry, const double& rz, const size_t& dir) const override{
        double res=0.;
        if(dir==0){
            auto const rx_ad=make_fvar<double,1>(rx);
            auto const autodiff=psi(rx_ad,ry,rz);
            res=autodiff.derivative(1);
        }
        else if(dir==1){
            auto const ry_ad=make_fvar<double,1>(ry);
            auto const autodiff=psi(rx,ry_ad,rz);
            res=autodiff.derivative(1);
        }
        else if(dir==2){
            auto const rz_ad=make_fvar<double,1>(rz);
            auto const autodiff=psi(rx,ry,rz_ad);
            res=autodiff.derivative(1);
        }
        return res;
    }
private:
    template<typename Tx, typename Ty, typename Tz>
    promote<Tx,Ty,Tz> psi(const Tx& rx, const Ty& ry, const Tz& rz) const{
        auto r=sqrt(rx*rx+ry*ry+rz*rz);
        return _norm*rz*exp(-0.5*_Z*r);//because cosT*r=rz
    }

    double _Z;
    double _norm;
};

/*class Pzorb_integr:
public Orbital{
public:
    Pzorb_integr(const double& Z):
    _Z(Z){}

    double operator()(const double& rx, const double& ry, const double& rz) const override{
        double r=sqrt(rx*rx+ry*ry+rz*rz);
        return -2.*exp(-0.5*_Z*r)*(2.+r*_Z)/_Z;
    }
private:
    double _Z;
};*/

class Atom{
    friend class AtomsSet;
    friend class Material;
public:
    Atom(const double& x, const double& y, const double& z,
         Orbital* orb):
    _x(x),_y(y),_z(z),
    _orb(orb){}
private:
    Orbital* _orb;
    double _x,_y,_z;
};

class AtomsSet{
public:
    AtomsSet(){}

    void add_atom(const Atom& atom){
        _atoms.push_back(atom);
    }

    Atom& get_atom(const size_t& ia){
        return _atoms[ia];
    }

    void compute_on_grid(Grid2D* kxygrid){
        size_t Nkx=kxygrid->size1();
        size_t Nky=kxygrid->size2();
        size_t Na=this->natoms();

        _indx=new MultiIndex({Nkx,Nky,Na});
        _exps=new complex_t[_indx->size()];

        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                for(size_t ia=0; ia<Na; ia++){
                    double kx=(*kxygrid)(ikx,iky)[0];
                    double ky=(*kxygrid)(ikx,iky)[1];

                    double Rx=_atoms[ia]._x;
                    double Ry=_atoms[ia]._y;

                    double k_R=kx*Rx+ky*Ry;

                    size_t indx=(*_indx)({ikx,iky,ia});
                    _exps[indx]=exp(I*k_R);
                }
            }
        }
        return;
    }

    complex_t Phi(const double& rx, const double& ry, const double& rz,
                  const double& kx, const double& ky) const{
        complex_t res=0.;
        for(auto atom: _atoms){
            double Rx=atom._x;
            double Ry=atom._y;
            double Rz=atom._z;

            double k_R=kx*Rx+ky*Ry;
            res+=exp(I*k_R)*atom._orb->p(rx-Rx,ry-Ry,rz-Rz);
        }
        return res;
    }

    complex_t Phi(const double& rx, const double& ry, const double& rz,
                  const size_t& ikx, const size_t& iky) const{
        complex_t res=0.;
        for(size_t ia=0; ia<this->natoms(); ia++){
            double Rx=_atoms[ia]._x;
            double Ry=_atoms[ia]._y;
            double Rz=_atoms[ia]._z;

            size_t indx=(*_indx)({ikx,iky,ia});
            res+=_exps[indx]*_atoms[ia]._orb->p(rx-Rx,ry-Ry,rz-Rz);
        }
        return res;
    }

    complex_t dPhi(const double& rx, const double& ry, const double& rz,
                   const size_t& ikx, const size_t& iky,
                   const size_t& dir) const{
        complex_t res=0.;
        for(size_t ia=0; ia<this->natoms(); ia++){
            double Rx=_atoms[ia]._x;
            double Ry=_atoms[ia]._y;
            double Rz=_atoms[ia]._z;

            size_t indx=(*_indx)({ikx,iky,ia});
            res+=_exps[indx]*_atoms[ia]._orb->dp(rx-Rx,ry-Ry,rz-Rz,dir);
        }
        return res;
    }

    size_t natoms() const{
        return _atoms.size();
    }
private:
    std::vector<Atom> _atoms;

    MultiIndex* _indx;
    complex_t* _exps;
};

class Material{
public:
    Material(){}
    ~Material(){}

    void add_atomsset(const AtomsSet& set){
        _sets.push_back(set);
    }

    template<typename T>
    complex_t PhiI(const size_t& i,
                   const double& rx, const double& ry, const double& rz,
                   const T& kx, const T& ky) const{
        return _sets[i].Phi(rx,ry,rz,kx,ky);
    }

    template<typename T>
    complex_t dPhiI(const size_t& i,
                    const double& rx, const double& ry, const double& rz,
                    const T& kx, const T& ky,
                    const size_t& dir) const{
        return _sets[i].dPhi(rx,ry,rz,kx,ky,dir);
    }

    void print_atoms(std::ofstream& out) const{
        out<<std::fixed;
        out<<std::setprecision(5);
        
        out<<this->natoms()<<std::endl;

        for(auto set: _sets){
            for(size_t ia=0; ia<set.natoms(); ia++){
                auto atom=set.get_atom(ia);
                
                out<<std::setw(12)<<atom._x;
                out<<std::setw(12)<<atom._y;
                out<<std::setw(12)<<atom._z;
                out<<std::endl;
            }
        }

        return;
    }

    size_t natoms() const{
        size_t n=0;
        for(auto set: _sets){
            n+=set.natoms();
        }
        return n;
    }
private:
    std::vector<AtomsSet> _sets;
};

#include <boost/numeric/ublas/matrix.hpp>

AtomsSet GenerateGraphenePattern(Orbital* orb,
    const double& a, const size_t& Nclx, const size_t& Ncly,
    const double& dx, const double& dy, const double& dz){
    
    typedef boost::numeric::ublas::vector<double> vector_t;

    AtomsSet set;
    
    vector_t a1(2);
    a1(0)=a/2.*sqrt(3.);
    a1(1)=a/2.;

    vector_t a2(2);
    a2(0)=a/2.*sqrt(3.);
    a2(1)=-a/2.;

    vector_t zero(2);
    zero(0)=dx;
    zero(1)=dy;

    for(int ix=0; ix<=2*Nclx; ix++){
        for(int iy=0; iy<=2*Ncly; iy++){
            int ixx=ix-Nclx;
            int iyy=iy-Ncly;

            vector_t posA1=zero+ixx*a1+iyy*a2;

            double A1_x=posA1[0];
            double A1_y=posA1[1];

            set.add_atom(Atom(A1_x,A1_y,dz,orb));
        }
    }

    return set;
}

#endif