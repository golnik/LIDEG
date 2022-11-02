#include <iostream>
#include <vector>

#include "mini/ini.h"

#include "utils/utils.hpp"
#include "utils/grid.hpp"
#include "utils/multiarray.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "model1/graphene.hpp"
#include "model1/WFs.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/format.hpp>

typedef matrix<double> matrix_t;

#define Nx_max 64
#define Ny_max 64
#define Nz_max 64

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];
        int tstep=std::stoi(argv[2])-1;

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

        ExternalField* E0=nullptr;

        //prepare xyz grids
        /*vector<double> a1(2);
        vector<double> a2(2);
        vector<double> origin(2);

        if(params.rgrid_type==regular){
            origin(0)=params.xmin;
            origin(1)=params.ymin;

            a1(0)=(params.xmax-params.xmin)/(params.Nx-1);
            a1(1)=0.;

            a2(0)=0.;
            a2(1)=(params.ymax-params.ymin)/(params.Ny-1);
        }
        else if(params.rgrid_type==ucell){
            origin(0)=-(1./sqrt(3.))*params.a;
            origin(1)=0.0;

            a1(0)=params.a/2.*sqrt(3.);
            a1(1)=params.a/2.;
            a1/=(params.Nx-1);
            
            a2(0)=params.a/2.*sqrt(3.);
            a2(1)=-params.a/2.;
            a2/=(params.Ny-1);
        }

        //write grids to file
        std::ofstream grid_out(params.rgfile_fname);
        grid_out<<params.Nx*params.Ny<<std::endl;

        grid_out<<std::scientific;
        grid_out<<std::setprecision(8);

        MultiArray<std::pair<double,double>,Nx_max,Ny_max> xygrid;
        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                auto o=origin+ix*a1+iy*a2;
                
                //double x=xgrid[ix];
                //double y=ygrid[iy];
                double x=o(0);
                double y=o(1);

                xygrid(ix,iy)=std::make_pair(x,y);
                grid_out<<std::setw(20)<<x;
                grid_out<<std::setw(20)<<y<<std::endl;

                //o+=a2;
            }
            //o+=a1;
        }
        grid_out.close();*/

        double Ox=-(1./sqrt(3.))*params.a;
        double Oy=0.;
        double a1x=params.a/2.*sqrt(3.);
        double a1y=params.a/2.;
        double a2x=a1x;
        double a2y=-a1y;

        Grid2D* xygrid=new UCellGrid2D(Ox,Oy,a1x,a1y,params.Nx,a2x,a2y,params.Ny);

        //write grids to file
        std::ofstream grid_out(params.rgfile_fname);
        grid_out<<params.Nx*params.Ny<<std::endl;

        grid_out<<std::scientific;
        grid_out<<std::setprecision(8);

        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                grid_out<<std::setw(20)<<(*xygrid)(ix,iy)[0];
                grid_out<<std::setw(20)<<(*xygrid)(ix,iy)[1]<<std::endl;
            }
        }
        grid_out.close();

        Grid1D* zgrid;
        Pzorb* pz;
        if(params.Nz==0 || params.Nz==1){
            zgrid=new RegularGrid1D(0.0,1.0,1);
            pz=new Pzorb_integr(params.Z);
        }
        else{
            zgrid=new RegularGrid1D(params.zmin,params.zmax,params.Nz);
            pz=new Pzorb_normal(params.Z);
        }
        
        //create graphene model
        GrapheneModel gm(params.a,params.e2p,params.gamma,params.s,params.Td,
                         E0,E0);

        //create graphene layer
        double R0x=0.5*(params.xmax+params.xmin);
        double R0y=0.5*(params.ymax+params.ymin);        
        GrapheneLayer gl(params.a,
                params.Nclx,params.Ncly,
                R0x,R0y,params.Rmax);

        //print atom positions to file
        std::ofstream atoms_out;
        atoms_out.open(params.afile_fname);
        gl.print_atoms(atoms_out);
        atoms_out.close();

        //area of BZ
        double SBZ=pow(2.*M_PI,2.)/(0.5*sqrt(3.)*params.a*params.a);

        //create Dirac points
        //std::vector<double> Dirac_Kx;
        //std::vector<double> Dirac_Ky;
        //std::vector<int> Dirac_type;

        //gm.get_Dirac_points(Dirac_Kx,Dirac_Ky,Dirac_type);
        //size_t nK=Dirac_Kx.size();

        //read density from file
        std::string dens_t_fname=params.densfile_fname;
        replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (tstep+1)));

        std::cout<<"Reciprocal space density file: "<<dens_t_fname<<std::endl;

        matrix_t dens_vv=read_2D_from_file<matrix_t>(dens_t_fname,0,params.Nkx,params.Nky);
        matrix_t dens_cc=read_2D_from_file<matrix_t>(dens_t_fname,1,params.Nkx,params.Nky);
        matrix_t dens_cv_re=read_2D_from_file<matrix_t>(dens_t_fname,2,params.Nkx,params.Nky);
        matrix_t dens_cv_im=read_2D_from_file<matrix_t>(dens_t_fname,3,params.Nkx,params.Nky);

        //array for real space data
        MultiArray<double,Nx_max,Ny_max,Nz_max> res;
        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<zgrid->size(); iz++){
                    res(ix,iy,iz)=0.;
                }
            }
        }

        //create kgrid
        double Okx=0.;
        double Oky=0.;
        double b1x=2.*M_PI/(sqrt(3.)*params.a);
        double b1y=2.*M_PI/(params.a);
        double b2x=b1x;
        double b2y=-b1y;

        /*double Okx=0.;
        double Oky=M_PI;
        double b1x=M_PI-Okx;
        double b1y=2.*M_PI-Oky;
        double b2x=M_PI-Okx;
        double b2y=0.-Oky;*/

        Grid2D* kxygrid=new UCellGrid2D(Okx,Oky,b1x,b1y,params.Nkx,b2x,b2y,params.Nky);

        WFs wfs(&gm,&gl,pz);
        WFs_grid wfs_g(&gm,&gl,pz,kxygrid);

        //create multi grid
        //size_t Nmulti=params.Nkx*params.Nky;
        //double* mgrid=new double [2*Nmulti];
        //double* mres=new double [Nmulti];

        //size_t indx=0;
        //for(size_t ikx=0; ikx<params.Nkx; ikx++){
        //    for(size_t iky=0; iky<params.Nky; iky++){
        //        mgrid[indx         ]=ikx;
        //        mgrid[indx+Nmulti  ]=iky;
        //        indx++;
        //    }
        //}

        /*double progress=0.;
        double dP=1./(static_cast<double>(params.Nx)*static_cast<double>(nK));
        for(size_t iK=0; iK<nK; iK++){//loop over Dirac points
            double Kx=Dirac_Kx[iK];
            double Ky=Dirac_Ky[iK];
            double Ktype=Dirac_type[iK];

            double kxmin=Kx-params.dkx;
            double kxmax=Kx+params.dkx;
            double kymin=Ky-params.dky;
            double kymax=Ky+params.dky;

            auto kx_grid=create_grid(kxmin,kxmax,params.Nkx,params.kgrid_type);
            auto ky_grid=create_grid(kymin,kymax,params.Nky,params.kgrid_type);

            //create WFs model
            WFs wfs(&gm,&gl,params.Z);
            WFs_grid wfs_g(&gm,&gl,params.Z,kx_grid,ky_grid);

            //double dkx=kx_grid[1]-kx_grid[0];
            //double dky=ky_grid[1]-ky_grid[0];

            //3D real-space densities
            for(size_t ix=0; ix<params.Nx; ix++){
                progress+=dP;
                if((ix%static_cast<int>(2*log2(params.Nx)))==0){
                    printProgress(progress);
                }

                for(size_t iy=0; iy<params.Ny; iy++){
                    for(size_t iz=0; iz<params.Nz; iz++){
                        //double x=xgrid[ix];
                        //double y=ygrid[iy];
                        //double x=xygrid(ix,iy).first;
                        //double y=xygrid(ix,iy).second;
                        double x=(*xygrid)(ix,iy)[0];
                        double y=(*xygrid)(ix,iy)[1];
                        double z=zgrid[iz];

                        res(ix,iy,iz)+=integrate(params.Nkx,params.Nky,
                            [params,
                            dens_vv,dens_cc,dens_cv_re,dens_cv_im,
                            Ktype,SBZ,
                            &wfs,&wfs_g,
                            kx_grid,ky_grid,
                            x,y,z](const size_t& ikx_K, const size_t& iky_K){
                                size_t ikx=0;
                                size_t iky=0;
                                if(Ktype==0){//normal indices for K point
                                    ikx=ikx_K;
                                    iky=iky_K;
                                }
                                else{//inverse indices for K' point
                                    ikx=ikx_K;
                                    iky=params.Nky-iky_K-1;
                                }

                                complex_t dens_cv=dens_cv_re(ikx,iky)+I*dens_cv_im(ikx,iky);

                                complex_t psip=wfs_g.psip(x,y,z,ikx_K,iky_K);
                                complex_t psim=wfs_g.psim(x,y,z,ikx_K,iky_K);

                                double rho_vv=pow(std::abs(psip),2.);
                                double rho_cc=pow(std::abs(psim),2.);
                                complex_t rho_vc=std::conj(psip)*psim;

                                double rho_t=dens_vv(ikx,iky)*rho_vv
                                            +dens_cc(ikx,iky)*rho_cc;
                                            //+2.*std::real(dens_cv*rho_vc);

                                double res=2.*3.*(2./SBZ)*(rho_t-rho_vv);

                                //double res=rho_cc-rho_vv;

                                //double res=std::abs(wfs_g.psip(x,y,z,ikx_K,iky_K))
                                //          -std::abs(wfs.psip(x,y,z,kx_grid[ikx],ky_grid[iky]));

                                return res;
                            },
                            kxmin,kxmax,
                            kymin,kymax                    
                        );
                    }//z loop
                }//y loop        
            }//x loop
        }//loop over Dirac points        
        */

        //3D real-space densities
        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<zgrid->size(); iz++){
                    double x=(*xygrid)(ix,iy)[0];
                    double y=(*xygrid)(ix,iy)[1];
                    double z=(*zgrid)[iz];

                    auto func=[x,y,z,
                    dens_vv,dens_cc,dens_cv_re,dens_cv_im,
                    &wfs,&wfs_g,
                    SBZ,
                    kxygrid](const size_t& ikx, const size_t& iky){
                        double kx=(*kxygrid)(ikx,iky)[0];
                        double ky=(*kxygrid)(ikx,iky)[1];

                        //complex_t psip=wfs.psip(x,y,z,kx,ky);
                        //complex_t psim=wfs.psim(x,y,z,kx,ky);

                        complex_t psip=wfs_g.psip(x,y,z,ikx,iky);
                        complex_t psim=wfs_g.psim(x,y,z,ikx,iky);

                        double rho_vv=pow(std::abs(psip),2.);
                        double rho_cc=pow(std::abs(psim),2.);
                        complex_t rho_vc=std::conj(psip)*psim;

                        complex_t dens_cv=dens_cv_re(ikx,iky)+I*dens_cv_im(ikx,iky);

                        double rho_t=dens_vv(ikx,iky)*rho_vv
                                    +dens_cc(ikx,iky)*rho_cc
                                    +2.*std::real(dens_cv*rho_vc);

                        double res=(2./SBZ)*(rho_t-rho_vv);

                        //double res=(2./SBZ)*std::real(dens_cv*rho_vc);

                        //double res=(2./SBZ)*rho_t;

                        //double res=std::real(psim)-std::real(psim_g);

                        //double res=std::norm(psip)-std::norm(psim);
                        //double res=std::norm(psip_g)-std::norm(psim_g);

                        //double res=kx*kx+ky*ky;

                        //double res=pow(kx-ky,2.)*pow(cos(kx+ky),2.);

                        return res;
                    };

                    //double integr=kxygrid->integrate(func);
                    //std::cout<<integr<<std::endl;

                    res(ix,iy,iz)+=kxygrid->integrate(func);
                }
            }
        }

        //output real space data
        std::string rho_t_fname=params.rhofile_fname;
        replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (tstep+1)));
        std::ofstream rho_t_out(rho_t_fname);

        std::cout<<"Real space density will be written to: "<<rho_t_fname<<std::endl;

        rho_t_out<<std::scientific;
        rho_t_out<<std::setprecision(10);
        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<zgrid->size(); iz++){
                    rho_t_out<<std::setw(20)<<res(ix,iy,iz)<<std::endl;
                }
            }
        }
        rho_t_out.close();

        double sum=0.;
        for(size_t iz=0; iz<zgrid->size(); iz++){
            std::function<double(const size_t& ix, const size_t& iy)> rho_xy=
            [res,iz](const size_t& ix, const size_t& iy){
                return res(ix,iy,iz);
            };
            sum+=xygrid->integrate(rho_xy);
        }
        std::cout<<"Sum: "<<sum<<std::endl;

    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
