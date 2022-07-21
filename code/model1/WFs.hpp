#ifndef WFS_HPP
#define WFS_HPP

#include "utils.hpp"
#include "graphene.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class WFs{
    typedef vector<double> vector_t;
public:
    WFs(GrapheneModel* gm, const double& Z,
        const size_t& Nclx, const size_t& Ncly):
    _gm{gm},_Z(Z){
        double a=gm->_a;

        vector_t a1(2);
        a1(0)=a/2.*sqrt(3.);
        a1(1)=a/2.;

        vector_t a2(2);
        a2(0)=a/2.*sqrt(3.);
        a2(1)=-a/2.;

        vector_t diffA1A2(2);
        diffA1A2(0)=-a/sqrt(3.);

        vector_t zero(2);

        std::vector<vector_t> _A1;
        std::vector<vector_t> _A2;

        for(int ix=0; ix<=2*Nclx; ix++){
            for(int iy=0; iy<=2*Ncly; iy++){
                int ixx=ix-Nclx;
                int iyy=iy-Ncly;

                vector_t posA1=zero+ixx*a1+iyy*a2;
                vector_t posA2=posA1-diffA1A2;

                _A1.push_back(posA1);
                _A2.push_back(posA2);
            }
        }

        _Natoms=_A1.size();
        _A1_x=new double [_Natoms];
        _A1_y=new double [_Natoms];
        _A2_x=new double [_Natoms];
        _A2_y=new double [_Natoms];

        for(size_t ia=0; ia<_Natoms; ia++){
            _A1_x[ia]=_A1[ia][0];
            _A1_y[ia]=_A1[ia][1];
            _A2_x[ia]=_A2[ia][0];
            _A2_y[ia]=_A2[ia][1];
        }

    }

    ~WFs(){}

    double phi_2pz(const double& rx, const double& ry, const double& rz) const{
        double r=sqrt(rx*rx+ry*ry+rz*rz);
        double cos_theta=rz/r;

        //return _Z*r*cos_theta*exp(-0.5*_Z*r);

        return phi_2pz_intz(rx,ry,rz);
    }

    double phi_2pz_intz(const double& rx, const double& ry, const double& rz) const{
        double r=sqrt(rx*rx+ry*ry+rz*rz);
        return -2.*exp(-0.5*_Z*r)*(2.+r*_Z)/_Z;
    }

    complex_t PhiA(const double& rx, const double& ry, const double& rz,
                   const double& kx, const double& ky,
                   double* Ax, double* Ay) const{
        complex_t res=0.;
        
        for(size_t ia=0; ia<_Natoms; ia++){
            double Rx=Ax[ia];
            double Ry=Ay[ia];

            double k_R=kx*Rx+ky*Ry;
            res+=exp(I*k_R)*this->phi_2pz(rx-Rx,ry-Ry,rz);
        }
        return res;
    }

    complex_t PhiA1(const double& rx, const double& ry, const double& rz,
                    const double& kx, const double& ky) const{
        return PhiA(rx,ry,rz,kx,ky,_A1_x,_A1_y);
    }

    complex_t PhiA2(const double& rx, const double& ry, const double& rz,
                    const double& kx, const double& ky) const{
        return PhiA(rx,ry,rz,kx,ky,_A2_x,_A2_y);
    }

    complex_t psip(const double& rx, const double& ry, const double& rz,
                    const double& kx, const double& ky) const{
        return 1./(sqrt(2.*(1.+_gm->_s*abs(_gm->f(kx,ky)))))
            *(PhiA1(rx,ry,rz,kx,ky)+exp(-I*_gm->phi(kx,ky))*PhiA2(rx,ry,rz,kx,ky));
    }

    complex_t psim(const double& rx, const double& ry, const double& rz,
                    const double& kx, const double& ky) const{
        return 1./(sqrt(2.*(1.-_gm->_s*abs(_gm->f(kx,ky)))))
            *(PhiA1(rx,ry,rz,kx,ky)-exp(-I*_gm->phi(kx,ky))*PhiA2(rx,ry,rz,kx,ky));
    }    
//private:
    GrapheneModel* _gm;
    double _Z;

    //std::vector<vector_t> _A1;
    //std::vector<vector_t> _A2;

    size_t _Natoms;
    double* _A1_x;
    double* _A1_y;
    double* _A2_x;
    double* _A2_y;
};

#endif