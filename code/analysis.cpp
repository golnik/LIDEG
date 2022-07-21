#include <iostream>
#include <vector>

#include "mini/ini.h"

#include "utils.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "model1/graphene.hpp"
#include "model1/WFs.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/format.hpp>

typedef matrix<double> matrix_t;

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

        //params.xmin=1./params.kx_min;
        //params.xmax=1./params.kx_max;
        //params.ymin=1./params.ky_min;
        //params.ymax=1./params.ky_max;

        //prepare grids
        auto tgrid  =create_grid(params.tmin,params.tmax,params.Nt);
        auto kx_grid=create_grid(params.kx_min,params.kx_max,params.Nkx);
        auto ky_grid=create_grid(params.ky_min,params.ky_max,params.Nky);

        //prepare xyz grids
        auto xgrid=create_grid(params.xmin,params.xmax,params.Nx);
        auto ygrid=create_grid(params.ymin,params.ymax,params.Ny);
        auto zgrid=create_grid(params.zmin,params.zmax,params.Nz);

        double dz=zgrid[1]-zgrid[0];
        double dkx=kx_grid[1]-kx_grid[0];
        double dky=ky_grid[1]-ky_grid[0];

        double z=params.zmax;

        //create multi grid
        size_t Nmulti=params.Nkx*params.Nky;
        double* mgrid=new double [2*Nmulti];

        size_t indx=0;
        for(size_t ikx=0; ikx<params.Nkx; ikx++){
            for(size_t iky=0; iky<params.Nky; iky++){
                mgrid[indx         ]=ikx;
                mgrid[indx+Nmulti  ]=iky;
                indx++;
            }
        }

        ExternalField* E0=nullptr;

        //create graphene model
        GrapheneModel gm(params.a,params.e2p,params.gamma,params.s,
                         E0,E0);

        //create WFs model
        WFs wfs(&gm,params.Z,params.Nclx,params.Ncly);

        //auto vec=wfs._A1;
        //std::cout<<vec.size()<<std::endl;

        //return 0;

        for(size_t it=0; it<params.Nt; it++){
            //read density from file
            std::string rho_t_fname=params.rhofile_fname;
            replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));

            matrix_t rho00=read_2D_from_file<matrix_t>(rho_t_fname,0,params.Nkx,params.Nky);
            matrix_t rho11=read_2D_from_file<matrix_t>(rho_t_fname,1,params.Nkx,params.Nky);
            matrix_t rho01=read_2D_from_file<matrix_t>(rho_t_fname,2,params.Nkx,params.Nky);

            //prepare output streams
            std::string dens_t_fname=params.densfile_fname;
            replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));
            std::ofstream dens_t_out(dens_t_fname);

            std::cout<<dens_t_fname<<std::endl;

            matrix<double> res(params.Nx,params.Ny);
            for(size_t ix=0; ix<params.Nx; ix++){
                for(size_t iy=0; iy<params.Ny; iy++){
                    res(ix,iy)=0.;
                }
            }

            for(size_t ix=0; ix<params.Nx; ix++){
                std::cout<<ix<<std::endl;     
                for(size_t iy=0; iy<params.Ny; iy++){
                    double x=xgrid[ix];
                    double y=ygrid[iy];
                    
                    double* mres=new double [Nmulti];

                    #pragma omp parallel for
                    for(size_t im=0; im<Nmulti; im++){
                        size_t ikx=mgrid[im];
                        size_t iky=mgrid[im+Nmulti];

                        double kx=kx_grid[ikx];
                        double ky=ky_grid[iky];

                        //mres[im]=abs(wfs.psip(x,y,z,kx,ky));
                        //mres[im]=std::real(std::conj(wfs.psim(x,y,z,kx,ky))*wfs.psip(x,y,z,kx,ky));

                        mres[im]=pow(std::abs(wfs.psip(x,y,z,kx,ky)),2.)*rho00(ikx,iky)
                                +pow(std::abs(wfs.psim(x,y,z,kx,ky)),2.)*rho11(ikx,iky)
                                +2.*std::real(std::conj(wfs.psim(x,y,z,kx,ky))*wfs.psip(x,y,z,kx,ky))*rho01(ikx,iky)
                                ;
                    }

                    for(size_t im=0; im<Nmulti; im++){
                        res(ix,iy)+=mres[im];
                    }

                    res(ix,iy)*=2.*dkx*dky;

                    /*std::vector<double> zvals(params.Nz);

                    #pragma omp parallel for
                    for(size_t iz=0; iz<params.Nz; iz++){

                        double x=xgrid[ix];
                        double y=ygrid[iy];
                        double z=zgrid[iz];

                        /*double int_kx=0.;
                        for(size_t ikx=0; ikx<params.Nkx; ikx++){
                            double kx=kx_grid[ikx];
                            
                            double int_ky=0.;
                            for(size_t iky=0; iky<params.Nky; iky++){
                                double ky=ky_grid[iky];

                                int_ky+=abs(wfs.PhiA1(x,y,z,kx,ky));
                            }
                            int_kx+=int_ky;
                        }

                        double int_kx=abs(wfs.psim(x,y,z,0,0));

                        zvals[iz]=int_kx;

                        //int_z+=pow(wfs.phi_2pz(x,y,z),2.);
                    }

                    double int_z=0.;
                    for(size_t iz=0; iz<params.Nz; iz++){
                        int_z+=zvals[iz];
                    }

                    res(ix,iy)=int_z;*/
                }
            }

            for(size_t ix=0; ix<params.Nx; ix++){
                double x=xgrid[ix];
                for(size_t iy=0; iy<params.Ny; iy++){
                    double y=ygrid[iy];



                    dens_t_out<<res(ix,iy)<<std::endl;
                }
            }

            dens_t_out.close();
        }
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
