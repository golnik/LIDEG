#ifndef WFS_HPP
#define WFS_HPP

#include <vector>

//#include "utils.hpp"
#include "graphene.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class GrapheneLayer{
    typedef vector<double> vector_t;
public:
    GrapheneLayer(const double& a,
        const size_t& Nclx, const size_t& Ncly,
        const double& R0x, const double& R0y,
        const double& Rmax){

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

                double A1_x=posA1[0];
                double A1_y=posA1[1];

                double R=sqrt(pow(A1_x-R0x,2.)+pow(A1_y-R0y,2.));
                if(R<=Rmax){
                    _A1.push_back(posA1);
                    _A2.push_back(posA2);
                }
            }
        }

        _Natoms=_A1.size();

        _A1_x=new double [_Natoms];
        _A1_y=new double [_Natoms];
        _A2_x=new double [_Natoms];
        _A2_y=new double [_Natoms];

        for(size_t ia=0; ia<_Natoms; ia++){
            _A1_x[ia]=_A1[ia](0);
            _A1_y[ia]=_A1[ia](1);
            _A2_x[ia]=_A2[ia](0);
            _A2_y[ia]=_A2[ia](1);
        }
    }

    ~GrapheneLayer(){
        delete[] _A1_x;
        delete[] _A1_y;
        delete[] _A2_x;
        delete[] _A2_y;
    }

    void print_atoms(std::ofstream& out) const{
        out<<std::fixed;
        out<<std::setprecision(5);
        
        out<<_Natoms<<std::endl;
        for(size_t ia=0; ia<_Natoms; ia++){
            out<<std::setw(12)<<_A1_x[ia];
            out<<std::setw(12)<<_A1_y[ia];
            out<<std::setw(12)<<_A2_x[ia];
            out<<std::setw(12)<<_A2_y[ia];
            out<<std::endl;
        }

        return;
    }

    size_t get_Natoms() const{
        return _Natoms;
    }

    void get_xy_A1(const size_t& ia, double& x, double& y) const{
        x=_A1_x[ia];
        y=_A1_y[ia];
        return;
    }

    void get_xy_A2(const size_t& ia, double& x, double& y) const{
        x=_A2_x[ia];
        y=_A2_y[ia];
        return;
    }    
private:
    size_t _Natoms;
    double* _A1_x;
    double* _A1_y;
    double* _A2_x;
    double* _A2_y;
};

class WFs{
public:
    WFs(GrapheneModel* gm, GrapheneLayer* gl,
        const double& Z):
    _gm{gm},_gl{gl},
    _Z(Z){}

    ~WFs(){
        //delete _gm;
        //delete _gl;
    }

    double phi_2pz(const double& rx, const double& ry, const double& rz) const{
        double r=sqrt(rx*rx+ry*ry+rz*rz);
        double cos_theta=rz/r;

        return _Z*r*cos_theta*exp(-0.5*_Z*r);

        //return phi_2pz_intz(rx,ry,rz);
    }

    double phi_2pz_intz(const double& rx, const double& ry, const double& rz) const{
        double r=sqrt(rx*rx+ry*ry+rz*rz);
        return -2.*exp(-0.5*_Z*r)*(2.+r*_Z)/_Z;
    }

    complex_t PhiA1(const double& rx, const double& ry, const double& rz,
                    const double& kx, const double& ky) const{
        complex_t res=0.;
        
        for(size_t ia=0; ia<_gl->get_Natoms(); ia++){
            double Rx=0.;
            double Ry=0.;
            _gl->get_xy_A1(ia,Rx,Ry);

            double k_R=kx*Rx+ky*Ry;
            res+=exp(I*k_R)*this->phi_2pz(rx-Rx,ry-Ry,rz);
        }

        return res;
    }

    complex_t PhiA2(const double& rx, const double& ry, const double& rz,
                    const double& kx, const double& ky) const{
        complex_t res=0.;
        
        for(size_t ia=0; ia<_gl->get_Natoms(); ia++){
            double Rx=0.;
            double Ry=0.;
            _gl->get_xy_A2(ia,Rx,Ry);

            double k_R=kx*Rx+ky*Ry;
            res+=exp(I*k_R)*this->phi_2pz(rx-Rx,ry-Ry,rz);
        }

        return res;
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
protected:
    GrapheneModel* _gm;
    GrapheneLayer* _gl;
    double _Z;
};

class WFs_grid:
public WFs{
public:
    WFs_grid(GrapheneModel* gm, GrapheneLayer* gl,
             const double& Z,
             const std::vector<double>& kx_grid, const std::vector<double>& ky_grid):
    WFs(gm,gl,Z),
    _kx_grid(kx_grid),_ky_grid(ky_grid){
        size_t Nkx=_kx_grid.size();
        size_t Nky=_ky_grid.size();

        size_t Natoms=WFs::_gl->get_Natoms();

        _exps_A1=new complex_t [Natoms*Nkx*Nky];
        _exps_A2=new complex_t [Natoms*Nkx*Nky];

        _pf_p_on_grid=new double [Nkx*Nky];
        _pf_m_on_grid=new double [Nkx*Nky];
        _exp_phi_on_grid=new complex_t [Nkx*Nky];

        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                double kx=kx_grid[ikx];
                double ky=ky_grid[iky];

                size_t ikxiky=indx_2D(ikx,iky);

                _pf_p_on_grid[ikxiky]=1./(sqrt(2.*(1.
                    +WFs::_gm->_s*abs(WFs::_gm->f(kx,ky)))));

                _pf_m_on_grid[ikxiky]=1./(sqrt(2.*(1.
                    -WFs::_gm->_s*abs(WFs::_gm->f(kx,ky)))));

                _exp_phi_on_grid[ikxiky]=exp(-I*_gm->phi(kx,ky));

                for(size_t ia=0; ia<Natoms; ia++){
                    size_t ikxikyia=indx_3D(ikx,iky,ia);

                    double Rx_A1=0.;
                    double Ry_A1=0.;
                    WFs::_gl->get_xy_A1(ia,Rx_A1,Ry_A1);

                    double k_R_A1=kx*Rx_A1+ky*Ry_A1;
                    _exps_A1[ikxikyia]=exp(I*k_R_A1);

                    double Rx_A2=0.;
                    double Ry_A2=0.;
                    WFs::_gl->get_xy_A2(ia,Rx_A2,Ry_A2);

                    double k_R_A2=kx*Rx_A2+ky*Ry_A2;
                    _exps_A2[ikxikyia]=exp(I*k_R_A2);
                }                  
            }
        }
    }

    ~WFs_grid(){
        delete[] _exps_A1;
        delete[] _exps_A2;
        delete[] _pf_p_on_grid;
        delete[] _pf_m_on_grid;
        delete[] _exp_phi_on_grid;
    }

    complex_t PhiA1(const double& rx, const double& ry, const double& rz,
                    const size_t& ikx, const size_t& iky) const{
        complex_t res=0.;
        
        for(size_t ia=0; ia<WFs::_gl->get_Natoms(); ia++){
            double Rx=0.;
            double Ry=0.;
            _gl->get_xy_A1(ia,Rx,Ry);

            size_t ikxikyia=indx_3D(ikx,iky,ia);
            res+=_exps_A1[ikxikyia]*this->phi_2pz(rx-Rx,ry-Ry,rz);
        }

        return res;
    }

    complex_t PhiA2(const double& rx, const double& ry, const double& rz,
                    const size_t& ikx, const size_t& iky) const{
        complex_t res=0.;
        
        for(size_t ia=0; ia<WFs::_gl->get_Natoms(); ia++){
            double Rx=0.;
            double Ry=0.;
            _gl->get_xy_A2(ia,Rx,Ry);

            size_t ikxikyia=indx_3D(ikx,iky,ia);
            res+=_exps_A2[ikxikyia]*this->phi_2pz(rx-Rx,ry-Ry,rz);
        }

        return res;
    }

    complex_t psip(const double& rx, const double& ry, const double& rz,
                   const size_t& ikx, const size_t& iky) const{
        size_t ikxiky=indx_2D(ikx,iky);
        return _pf_p_on_grid[ikxiky]
            *(PhiA1(rx,ry,rz,ikx,iky)+_exp_phi_on_grid[ikxiky]*PhiA2(rx,ry,rz,ikx,iky));
    }

    complex_t psim(const double& rx, const double& ry, const double& rz,
                   const size_t& ikx, const size_t& iky) const{
        size_t ikxiky=indx_2D(ikx,iky);                    
        return _pf_m_on_grid[ikxiky]
            *(PhiA1(rx,ry,rz,ikx,iky)-_exp_phi_on_grid[ikxiky]*PhiA2(rx,ry,rz,ikx,iky));
    }   
private:
    complex_t* _exps_A1;
    complex_t* _exps_A2;
    double* _pf_p_on_grid;
    double* _pf_m_on_grid;
    complex_t* _exp_phi_on_grid;

    std::vector<double> _kx_grid;
    std::vector<double> _ky_grid;

    size_t indx_2D(const size_t& i, const size_t& j) const{
        return i*_kx_grid.size()+j;
    }

    size_t indx_3D(const size_t& i, const size_t& j, const size_t& k) const{
        return k*_kx_grid.size()*_ky_grid.size()+indx_2D(i,j);
    }
};

#endif