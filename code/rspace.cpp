#include <iostream>
#include <vector>

#include "mini/ini.h"

#include "utils/utils.hpp"
#include "utils/multiarray.hpp"
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
        int tstep=std::stoi(argv[2])-1;

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

        ExternalField* E0=nullptr;

        //prepare xyz grids
        auto xgrid=create_grid(params.xmin,params.xmax,params.Nx);
        auto ygrid=create_grid(params.ymin,params.ymax,params.Ny);
        auto zgrid=create_grid(params.zmin,params.zmax,params.Nz);

        //create graphene model
        GrapheneModel gm(params.a,params.e2p,params.gamma,params.s,
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
        std::vector<double> Dirac_Kx;
        std::vector<double> Dirac_Ky;
        std::vector<int> Dirac_type;

        gm.get_Dirac_points(Dirac_Kx,Dirac_Ky,Dirac_type);
        size_t nK=Dirac_Kx.size();

        //read density from file
        std::string dens_t_fname=params.densfile_fname;
        replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (tstep+1)));

        std::cout<<"Reciprocal space density file: "<<dens_t_fname<<std::endl;

        matrix_t dens_vv=read_2D_from_file<matrix_t>(dens_t_fname,0,params.Nkx,params.Nky);
        matrix_t dens_cc=read_2D_from_file<matrix_t>(dens_t_fname,1,params.Nkx,params.Nky);
        matrix_t dens_cv_re=read_2D_from_file<matrix_t>(dens_t_fname,2,params.Nkx,params.Nky);
        matrix_t dens_cv_im=read_2D_from_file<matrix_t>(dens_t_fname,3,params.Nkx,params.Nky);

        //array for real space data
        MultiArray<double,64,64,64> res;
        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<params.Nz; iz++){
                    res(ix,iy,iz)=0.;
                }
            }
        }

        //create multi grid
        size_t Nmulti=params.Nkx*params.Nky;
        double* mgrid=new double [2*Nmulti];
        double* mres=new double [Nmulti];

        size_t indx=0;
        for(size_t ikx=0; ikx<params.Nkx; ikx++){
            for(size_t iky=0; iky<params.Nky; iky++){
                mgrid[indx         ]=ikx;
                mgrid[indx+Nmulti  ]=iky;
                indx++;
            }
        }

        double progress=0.;
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
            //WFs wfs(&gm,&gl,params.Z);
            WFs_grid wfs_g(&gm,&gl,params.Z,kx_grid,ky_grid);

            double dkx=kx_grid[1]-kx_grid[0];
            double dky=ky_grid[1]-ky_grid[0];

            //3D real-space densities
            for(size_t ix=0; ix<params.Nx; ix++){
                progress+=dP;
                if((ix%static_cast<int>(2*log2(params.Nx)))==0){
                    printProgress(progress);
                }

                for(size_t iy=0; iy<params.Ny; iy++){
                    for(size_t iz=0; iz<params.Nz; iz++){
                        double x=xgrid[ix];
                        double y=ygrid[iy];
                        double z=zgrid[iz];

                        for(size_t im=0; im<Nmulti; im++)
                            mres[im]=0.;

                        #pragma omp parallel for
                        for(size_t im=0; im<Nmulti; im++){

                            size_t ikx_K=mgrid[im];
                            size_t iky_K=mgrid[im+Nmulti];

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

                            //double kx=kx_grid[ikx_K];
                            //double ky=ky_grid[iky_K];

                            //double kxt=kx+(*Afield_x)(time);
                            //double kyt=ky+(*Afield_y)(time);

                            complex_t dens_cv=dens_cv_re(ikx,iky)+I*dens_cv_im(ikx,iky);

                            //mres[im]=wfs.phi_2pz(x,y,z);

                            //mres[im]=wfs.PhiA1(x,y,z,kx,ky);
                            //mres[im]=wfs.PhiA2(x,y,z,kx,ky);

                            //double rho_vv=std::norm(wfs.psip(x,y,z,kx,ky));
                            //double rho_cc=std::norm(wfs.psim(x,y,z,kx,ky));

                            //complex_t psip=wfs.psip(x,y,z,kxt,kyt);
                            //complex_t psim=wfs.psim(x,y,z,kxt,kyt);

                            //complex_t rho_vc=std::conj(wfs.psip(x,y,z,kx,ky))*wfs.psim(x,y,z,kx,ky);

                            complex_t psip=wfs_g.psip(x,y,z,ikx_K,iky_K);
                            complex_t psim=wfs_g.psim(x,y,z,ikx_K,iky_K);

                            double rho_vv=std::norm(psip);
                            double rho_cc=std::norm(psim);
                            complex_t rho_vc=std::conj(psip)*psim;

                            double rho_t=dens_vv(ikx,iky)*rho_vv
                                        +dens_cc(ikx,iky)*rho_cc
                                        +2.*std::real(dens_cv*rho_vc);
                            //rho_t*=(2./SBZ);

                            mres[im]=rho_t-rho_vv;

                            //mres[im]=std::norm(psip);
                            //mres[im]=std::norm(psim);

                            //mres[im]=dens_vv(ikx,iky)*rho_vv
                            //        +dens_cc(ikx,iky)*rho_cc
                            //        -rho_vv;
                            //mres[im]=dens_cc(ikx,iky)
                            //        *(-rho_vv+rho_cc);
                            
                            //mres[im]=std::norm(wfs.psim(x,y,z,kx,ky));
                            //mres[im]=std::real(std::conj(wfs.psim(x,y,z,kx,ky))*wfs.psip(x,y,z,kx,ky));

                            //mres[im]=pow(std::abs(wfs.psip(x,y,z,kx,ky)),2.)*rho00(ikx,iky)
                            //        +pow(std::abs(wfs.psim(x,y,z,kx,ky)),2.)*rho11(ikx,iky)
                            //        +2.*std::real(std::conj(wfs.psim(x,y,z,kx,ky))*wfs.psip(x,y,z,kx,ky))*rho01(ikx,iky)
                            //        ;
                        }

                        for(size_t im=0; im<Nmulti; im++){
                            res(ix,iy,iz)+=mres[im];
                        }
                    }//z loop
                }//y loop        
            }//x loop
        }//loop over Dirac points        

        //output real space data
        std::string rho_t_fname=params.rhofile_fname;
        replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (tstep+1)));
        std::ofstream rho_t_out(rho_t_fname);

        std::cout<<"Real space density will be written to: "<<rho_t_fname<<std::endl;

        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                double res_xy=0.;
                for(size_t iz=0; iz<params.Nz; iz++){
                    //double x=xgrid[ix];
                    //double y=ygrid[iy];
                    //double z=zgrid[iz];

                    //res_xy+=res(ix,iy,iz);
                    rho_t_out<<res(ix,iy,iz)<<std::endl;
                }
                //rho_t_out<<res_xy<<std::endl;
            }
        }
        rho_t_out.close();

    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
