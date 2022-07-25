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

        //prepare grids
        auto tgrid  =create_grid(params.tmin,params.tmax,params.Nt);

        //prepare xyz grids
        auto xgrid=create_grid(params.xmin,params.xmax,params.Nx);
        auto ygrid=create_grid(params.ymin,params.ymax,params.Ny);
        auto zgrid=create_grid(params.zmin,params.zmax,params.Nz);

        double dz=zgrid[1]-zgrid[0];
        double dkx=0.;
        double dky=0.;

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

        //read fields from file
        std::vector<double> tgrid_fit;
        std::vector<double> Adata_x_fit;
        std::vector<double> Adata_y_fit;
        std::vector<double> Edata_x_fit;
        std::vector<double> Edata_y_fit;

        read_column_from_file(params.field_fname,0,tgrid_fit);
        read_column_from_file(params.field_fname,1,Adata_x_fit);
        read_column_from_file(params.field_fname,2,Adata_y_fit);
        read_column_from_file(params.field_fname,3,Edata_x_fit);
        read_column_from_file(params.field_fname,4,Edata_y_fit);

        double t0_fit=tgrid_fit[0]/au2fs;
        double dt_fit=(tgrid_fit[1]-tgrid_fit[0])/au2fs;

        ExternalField* Afield_x=new ExternalFieldFromData(Adata_x_fit,t0_fit,dt_fit);
        ExternalField* Afield_y=new ExternalFieldFromData(Adata_y_fit,t0_fit,dt_fit);
        ExternalField* Efield_x=new ExternalFieldFromData(Edata_x_fit,t0_fit,dt_fit);
        ExternalField* Efield_y=new ExternalFieldFromData(Edata_y_fit,t0_fit,dt_fit);

        //create graphene model
        GrapheneModel gm(params.a,params.e2p,params.gamma,params.s,
                         Efield_x,Efield_y);

        //create WFs model
        WFs wfs(&gm,params.Z,params.Nclx,params.Ncly);

        //create Dirac points
        std::vector<double> Dirac_Kx;
        std::vector<double> Dirac_Ky;
        std::vector<int> Dirac_type;

        gm.get_Dirac_points(Dirac_Kx,Dirac_Ky,Dirac_type);
        size_t nK=Dirac_Kx.size();

        //auto vec=wfs._A1;
        //std::cout<<vec.size()<<std::endl;

        //return 0;

        //prepare output streams
        std::ofstream tfile_out;
        tfile_out.open(params.tfile_fname);

        tfile_out<<std::fixed;
        tfile_out<<std::setprecision(8);

        tfile_out<<"#";
        tfile_out<<std::setw(14)<<"Time";
        tfile_out<<std::setw(15)<<"Ax";
        tfile_out<<std::setw(15)<<"Ay";
        tfile_out<<std::setw(15)<<"Ex";
        tfile_out<<std::setw(15)<<"Ey";
        tfile_out<<std::setw(15)<<"rho_vv";
        tfile_out<<std::setw(15)<<"rho_cc";
        tfile_out<<std::setw(15)<<"Re{rho_cv}";
        tfile_out<<std::setw(15)<<"Im{rho_cv}";
        tfile_out<<std::setw(15)<<"Jx_intra";
        tfile_out<<std::setw(15)<<"Jy_intra";
        tfile_out<<std::setw(15)<<"Jx_inter";
        tfile_out<<std::setw(15)<<"Jy_inter";
        tfile_out<<std::endl;

        for(size_t it=0; it<params.Nt; it++){
            //read density from file
            std::string rho_t_fname=params.rhofile_fname;
            replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));

            matrix_t rho_vv=read_2D_from_file<matrix_t>(rho_t_fname,0,params.Nkx,params.Nky);
            matrix_t rho_cc=read_2D_from_file<matrix_t>(rho_t_fname,1,params.Nkx,params.Nky);
            matrix_t rho_cv_re=read_2D_from_file<matrix_t>(rho_t_fname,2,params.Nkx,params.Nky);
            matrix_t rho_cv_im=read_2D_from_file<matrix_t>(rho_t_fname,3,params.Nkx,params.Nky);

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

            double time=tgrid[it];

            //integrate over kgrids
            double pop0=0.;
            double pop1=0.;
            complex_t coh=0.;
            double Jra_x=0.;
            double Jra_y=0.;
            double Jer_x=0.;
            double Jer_y=0.;

            for(size_t iK=0; iK<nK; iK++){//loop over Dirac points
                double Kx=Dirac_Kx[iK];
                double Ky=Dirac_Ky[iK];

                auto kx_grid=create_grid(Kx-params.dkx,Kx+params.dkx,params.Nkx);
                auto ky_grid=create_grid(Ky-params.dky,Ky+params.dky,params.Nky);

                dkx=kx_grid[1]-kx_grid[0];
                dky=ky_grid[1]-ky_grid[0];

                for(size_t ikx_K=0; ikx_K<params.Nkx; ikx_K++){
                    for(size_t iky_K=0; iky_K<params.Nky; iky_K++){

                        size_t ikx=0;
                        size_t iky=0;
                        //normal indices for K point
                        if(Dirac_type[iK]==0){
                            ikx=ikx_K;
                            iky=iky_K;
                        }
                        else{
                            ikx=ikx_K;
                            iky=params.Nky-iky_K-1;
                        }

                        //integrate band populations
                        pop0+=rho_vv(ikx,iky);
                        pop1+=rho_cc(ikx,iky);

                        //integrate band coherences
                        complex_t rho_cv=rho_cv_re(ikx,iky)+I*rho_cv_im(ikx,iky);
                        coh+=rho_cv;

                        //integrate currents
                        double kxt=kx_grid[ikx_K]+(*Afield_x)(time);
                        double kyt=ky_grid[iky_K]+(*Afield_y)(time);

                        Jra_x+=gm.px_vv(kxt,kyt)*rho_vv(ikx,iky);
                              +gm.px_cc(kxt,kyt)*rho_cc(ikx,iky);

                        Jra_y+=gm.py_vv(kxt,kyt)*rho_vv(ikx,iky);
                              +gm.py_cc(kxt,kyt)*rho_cc(ikx,iky);

                        Jer_x+=2.*std::real(gm.px_cv(kxt,kyt)*rho_cv);
                        Jer_y+=2.*std::real(gm.py_cv(kxt,kyt)*rho_cv);
                    }
                }

                /*for(size_t ix=0; ix<params.Nx; ix++){
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
                    //}
                //}
            }//loop over Dirac points

            pop0 /= params.Nkx*params.Nky*nK;
            pop1 /= params.Nkx*params.Nky*nK;
            coh  /= params.Nkx*params.Nky*nK;

            Jra_x*=dkx*dky;
            Jra_y*=dkx*dky;
            Jer_x*=dkx*dky;
            Jer_y*=dkx*dky;

            tfile_out<<std::setw(15)<<time*au2fs;
            tfile_out<<std::setw(15)<<(*Afield_x)(time);
            tfile_out<<std::setw(15)<<(*Afield_y)(time);
            tfile_out<<std::setw(15)<<(*Efield_x)(time);
            tfile_out<<std::setw(15)<<(*Efield_y)(time);
            tfile_out<<std::setw(15)<<pop0;
            tfile_out<<std::setw(15)<<pop1;
            tfile_out<<std::setw(15)<<std::real(coh);
            tfile_out<<std::setw(15)<<std::imag(coh);
            tfile_out<<std::setw(15)<<Jra_x;
            tfile_out<<std::setw(15)<<Jra_y;
            tfile_out<<std::setw(15)<<Jer_x;
            tfile_out<<std::setw(15)<<Jer_y;
            tfile_out<<std::endl;

            for(size_t ix=0; ix<params.Nx; ix++){
                double x=xgrid[ix];
                for(size_t iy=0; iy<params.Ny; iy++){
                    double y=ygrid[iy];



                    dens_t_out<<res(ix,iy)<<std::endl;
                }
            }

            dens_t_out.close();
        }
        tfile_out.close();
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
