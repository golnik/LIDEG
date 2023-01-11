#ifndef GRID_HPP
#define GRID_HPP

#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/geometry/arithmetic/cross_product.hpp>

class Grid1D{
public:
    virtual size_t size() const=0;
    virtual double operator[](const size_t& i) const=0;
    virtual double get_dx() const=0;
private:
};

class RegularGrid1D:
public Grid1D{
public:
    RegularGrid1D(const double& min, const double& max, const size_t& N):
    _N(N){
        _grid=new double [N];

        if(N==1){
            _grid[0]=min;
            _dx=1.0;
        }
        else if(N>1){
            _dx=(max-min)/(N-1);
            double var=min;
            for(size_t i=0; i<N; i++){
                _grid[i]=var;
                var+=_dx;
            }
        }
    }

    virtual ~RegularGrid1D(){
        delete[] _grid;
    }

    size_t size() const override{
        return _N;
    }

    double operator[](const size_t& i) const override{
        return _grid[i];
    }

    double get_dx() const override{
        return _dx;
    }
private:
    double* _grid;
    double _dx;
    size_t _N;
};

/*class UniformGrid1D:
public Grid1D{
    typedef vector<double> vector_t;
public:    
    UniformGrid1D(const double& xmin, const double& xmax, const size_t& N,
                  const double& ymin=0., const double& ymax=0.):
    _N(N){
        _d=new vector_t(2);
        _O=new vector_t(2);

        (*_O)(0)=xmin;
        (*_O)(1)=ymin;

        (*_d)(0)=(xmax-xmin)/(N-1);
        (*_d)(1)=(ymax-ymin)/(N-1);
    }

    virtual ~UniformGrid1D(){
        delete _d;
        delete _O;
    }

    size_t size() const override{
        return _N;
    }

    std::array<double,2> operator[](const size_t& i) const override{
        auto res=(*_O)+i*(*_d);
        return {res[0],res[1]};
    }
private:
    vector_t* _d;
    vector_t* _O;
    size_t _N;
};*/


class GaussLegendreGrid1D:
public Grid1D{
public:
private:
};

class Grid2D{
public:
    typedef std::function<double(const size_t&, const size_t&)> func_t;

    virtual size_t size() const=0;
    virtual size_t size1() const=0;
    virtual size_t size2() const=0;

    virtual std::array<double,2> operator()(const size_t& i, const size_t& j) const=0;
    virtual double integrate(const func_t& f) const=0;

    virtual double get_dS() const=0;
private:
};

class RegularGrid2D:
public Grid2D{
public:
    RegularGrid2D(Grid1D* xgrid, Grid1D* ygrid):
    _xgrid{xgrid},_ygrid{ygrid}{}

    virtual ~RegularGrid2D(){
        delete _xgrid;
        delete _ygrid;
    }

    size_t size() const override{
        return _xgrid->size()*_ygrid->size();
    }

    size_t size1() const override{
        return _xgrid->size();
    }

    size_t size2() const override{
        return _ygrid->size();
    }

    std::array<double,2> operator()(const size_t& i, const size_t& j) const override{
        return {(*_xgrid)[i],(*_ygrid)[j]};
    }

    double integrate(const func_t& f) const override{
        return 0.;
    }

    double get_dS() const override{
        double dx=(*_xgrid)[1]-(*_xgrid)[0];
        double dy=(*_ygrid)[1]-(*_ygrid)[0];
        return dx*dy;
    }    
private:
    Grid1D* _xgrid;
    Grid1D* _ygrid;
};

class UCellGrid2D:
public Grid2D{
    typedef vector<double> vector_t;
public:
    UCellGrid2D(const double& Ox, const double& Oy,
                const double& x1, const double& y1, const size_t& Nx,
                const double& x2, const double& y2, const size_t& Ny):
    _Nx(Nx),_Ny(Ny){
        _a=new vector_t(2);
        _b=new vector_t(2);
        _O=new vector_t(2);

        (*_O)(0)=Ox;
        (*_O)(1)=Oy;

        (*_a)(0)=x1;
        (*_a)(1)=y1;
        (*_a)/=(Nx-1);

        (*_b)(0)=x2;
        (*_b)(1)=y2;
        (*_b)/=(Ny-1);
    }

    virtual ~UCellGrid2D(){
        delete _a;
        delete _b;
        delete _O;
    }

    size_t size() const override{
        return _Nx*_Ny;
    }

    size_t size1() const override{
        return _Nx;
    }

    size_t size2() const override{
        return _Ny;
    }    

    std::array<double,2> operator()(const size_t& i, const size_t& j) const override{
        auto res=(*_O)+i*(*_a)+j*(*_b);
        return {res(0),res(1)};
    }

    //2D trapezoidal integration
    double integrate(const func_t& f) const override{
        double res=0.;

        size_t Ncross=(_Nx-2)*(_Ny-2);
        double* tmp=new double[(_Nx-2)*(_Ny-2)];
        
        #pragma omp parallel for collapse(2)
        for(size_t ix=0; ix<_Nx-2; ix++){
            for(size_t iy=0; iy<_Ny-2; iy++){
                size_t indx=ix*(_Ny-2)+iy;

                tmp[indx]=f(ix+1,iy+1);
            }
        }
        for(size_t indx=0; indx<Ncross; indx++){
            res+=tmp[indx];
        }

        delete [] tmp;

        double sumi=0.;
        for(size_t ix=1; ix<_Nx-1; ix++){
            sumi+=f(ix,0)+f(ix,_Ny-1);
        }

        double sumj=0.;
        for(size_t iy=1; iy<_Ny-1; iy++){
            sumj+=f(0,iy)+f(_Nx-1,iy);
        }

        res+=0.5*(sumi+sumj)+0.25*(f(0,0)+f(_Nx-1,0)+f(0,_Ny-1)+f(_Nx-1,_Ny-1));

        return this->get_dS()*res;
    }

    double get_dS() const override{
        double ab=inner_prod(*_a,*_b);
        double na=norm_2(*_a);
        double nb=norm_2(*_b);
        double theta=acos(ab/(na*nb));
        double dS=na*na*sin(theta);

        return dS;
    }
private:
    size_t _Nx;
    size_t _Ny;

    vector_t* _a;
    vector_t* _b;
    vector_t* _O;    
};

class Integrator1D{
public:
    Integrator1D(Grid1D* grid):
    _grid(grid){}

    template<typename value_t, typename func_t>
    void trapz(const func_t& f, value_t& res) const{
        size_t N=_grid->size();
        double dx=_grid->get_dx();

        res=0.5*(f(0)+f(N-1));
        for(size_t i=1; i<N-1; i++){
            res+=f(i);
        }
        res*=dx;
        return;
    }
private:
    Grid1D* _grid;
};

class Integrator2D{
public:
    Integrator2D(Grid2D* grid):
    _grid{grid}{}

    template<typename value_t, typename func_t>
    void trapz(const func_t& f, value_t& res) const{
        size_t Nx=_grid->size1();
        size_t Ny=_grid->size2();
        double dS=_grid->get_dS();

        size_t Ncross=(Nx-2)*(Ny-2);
        value_t* tmp=new value_t[(Nx-2)*(Ny-2)];
        
        #pragma omp parallel for collapse(2)
        for(size_t ix=0; ix<Nx-2; ix++){
            for(size_t iy=0; iy<Ny-2; iy++){
                size_t indx=ix*(Ny-2)+iy;

                tmp[indx]=f(ix+1,iy+1);
            }
        }
        for(size_t indx=0; indx<Ncross; indx++){
            res+=tmp[indx];
        }

        delete [] tmp;

        for(size_t ix=1; ix<Nx-1; ix++){
            res+=0.5*(f(ix,0)+f(ix,Ny-1));
        }

        for(size_t iy=1; iy<Ny-1; iy++){
            res+=0.5*(f(0,iy)+f(Nx-1,iy));
        }

        res+=0.25*(f(0,0)+f(Nx-1,0)+f(0,Ny-1)+f(Nx-1,Ny-1));
        res*=dS;

        return;
    }
private:
    Grid2D* _grid;
};

#endif