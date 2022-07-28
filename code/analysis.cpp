#include <iostream>
#include <vector>

#include "mini/ini.h"

#include "utils/utils.hpp"
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

        ExternalField* Afield_x=new ExternalFieldFromData(Adata_x_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Afield_y=new ExternalFieldFromData(Adata_y_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Efield_x=new ExternalFieldFromData(Edata_x_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Efield_y=new ExternalFieldFromData(Edata_y_fit,t0_fit,dt_fit,params.E0);

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

        tfile_out<<std::scientific;
        //tfile_out<<std::fixed;
        tfile_out<<std::setprecision(8);

        tfile_out<<"#";
        tfile_out<<std::setw(19)<<"Time";
        tfile_out<<std::setw(20)<<"Ax";
        tfile_out<<std::setw(20)<<"Ay";
        tfile_out<<std::setw(20)<<"Ex";
        tfile_out<<std::setw(20)<<"Ey";
        tfile_out<<std::setw(20)<<"rho_vv";
        tfile_out<<std::setw(20)<<"rho_cc";
        tfile_out<<std::setw(20)<<"Re{rho_cv}";
        tfile_out<<std::setw(20)<<"Im{rho_cv}";
        tfile_out<<std::setw(20)<<"Jx_intra";
        tfile_out<<std::setw(20)<<"Jy_intra";
        tfile_out<<std::setw(20)<<"Jx_inter";
        tfile_out<<std::setw(20)<<"Jy_inter";
        tfile_out<<std::endl;

        for(size_t it=0; it<params.Nt; it++){
            //read density from file
            std::string dens_t_fname=params.densfile_fname;
            replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));

            std::cout<<dens_t_fname<<std::endl;

            matrix_t rho_vv=read_2D_from_file<matrix_t>(dens_t_fname,0,params.Nkx,params.Nky);
            matrix_t rho_cc=read_2D_from_file<matrix_t>(dens_t_fname,1,params.Nkx,params.Nky);
            matrix_t rho_cv_re=read_2D_from_file<matrix_t>(dens_t_fname,2,params.Nkx,params.Nky);
            matrix_t rho_cv_im=read_2D_from_file<matrix_t>(dens_t_fname,3,params.Nkx,params.Nky);

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
            double Jra[2]={0.,0.};
            double Jer[2]={0.,0.};

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

                dkx=kx_grid[1]-kx_grid[0];
                dky=ky_grid[1]-ky_grid[0];

                /*for(size_t ikx_K=0; ikx_K<params.Nkx; ikx_K++){
                    for(size_t iky_K=0; iky_K<params.Nky; iky_K++){

                        size_t ikx=0;
                        size_t iky=0;
                        //normal indices for K point
                        if(Ktype==0){
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

                        Jra[0]+=gm.px_vv(kxt,kyt)*rho_vv(ikx,iky);
                              +gm.px_cc(kxt,kyt)*rho_cc(ikx,iky);

                        Jra[1]+=gm.py_vv(kxt,kyt)*rho_vv(ikx,iky);
                              +gm.py_cc(kxt,kyt)*rho_cc(ikx,iky);

                        Jer[0]+=2.*std::real(gm.px_cv(kxt,kyt)*rho_cv);
                        Jer[1]+=2.*std::real(gm.py_cv(kxt,kyt)*rho_cv);
                    }
                }*/

                //integrate band populations
                pop0+=integrate(params.Nkx,params.Nky,
                    [rho_vv](const size_t& ix, const size_t& iy){
                        return rho_vv(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax
                );

                pop1+=integrate(params.Nkx,params.Nky,
                    [rho_cc](const size_t& ix, const size_t& iy){
                        return rho_cc(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax                    
                );                

                //integrate coherences
                double coh_re=integrate(params.Nkx,params.Nky,
                    [rho_cv_re](const size_t& ix, const size_t& iy){
                        return rho_cv_re(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax
                );

                double coh_im=integrate(params.Nkx,params.Nky,
                    [rho_cv_im](const size_t& ix, const size_t& iy){
                        return rho_cv_im(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax
                );

                coh+=coh_re+I*coh_im;

                for(int dir=0; dir<2; dir++){
                    Jra[dir]+=integrate(params.Nkx,params.Nky,
                        [rho_vv,rho_cc,
                        kx_grid,ky_grid,
                        time,
                        Afield_x,Afield_y,
                        Ktype,params,
                        gm,
                        dir](const size_t& ikx_K, const size_t& iky_K){
                            double kxt=kx_grid[ikx_K]+(*Afield_x)(time);
                            double kyt=ky_grid[iky_K]+(*Afield_y)(time);

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

                            if(dir==0){
                                return gm.px_vv(kxt,kyt)*rho_vv(ikx,iky)
                                      +gm.px_cc(kxt,kyt)*rho_cc(ikx,iky);
                            }
                            else if(dir==1){
                                return gm.py_vv(kxt,kyt)*rho_vv(ikx,iky)
                                      +gm.py_cc(kxt,kyt)*rho_cc(ikx,iky);
                            }
                        },
                        kxmin,kxmax,
                        kymin,kymax
                    );
                }

                for(int dir=0; dir<2; dir++){
                    Jer[dir]+=integrate(params.Nkx,params.Nky,
                        [rho_cv_re,rho_cv_im,
                        kx_grid,ky_grid,
                        time,
                        Afield_x,Afield_y,
                        Ktype,params,
                        gm,
                        dir](const size_t& ikx_K, const size_t& iky_K){
                            double kxt=kx_grid[ikx_K]+(*Afield_x)(time);
                            double kyt=ky_grid[iky_K]+(*Afield_y)(time);

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

                            complex_t rho_cv=rho_cv_re(ikx,iky)+I*rho_cv_im(ikx,iky);

                            if(dir==0){
                                return 2.*std::real(gm.px_cv(kxt,kyt)*rho_cv);
                            }
                            else if(dir==1){
                                return 2.*std::real(gm.py_cv(kxt,kyt)*rho_cv);
                            }
                        },
                        kxmin,kxmax,
                        kymin,kymax
                    );
                }                

                //res*=0.5*(b-a)*0.5*(d-c);

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

            /*pop0 /= params.Nkx*params.Nky*nK;
            pop1 /= params.Nkx*params.Nky*nK;
            coh  /= params.Nkx*params.Nky*nK;

            Jra[0]*=dkx*dky;
            Jra[1]*=dkx*dky;
            Jer[0]*=dkx*dky;
            Jer[1]*=dkx*dky;*/

            tfile_out<<std::setw(20)<<time*au2fs;
            tfile_out<<std::setw(20)<<(*Afield_x)(time);
            tfile_out<<std::setw(20)<<(*Afield_y)(time);
            tfile_out<<std::setw(20)<<(*Efield_x)(time);
            tfile_out<<std::setw(20)<<(*Efield_y)(time);
            tfile_out<<std::setw(20)<<pop0;
            tfile_out<<std::setw(20)<<pop1;
            tfile_out<<std::setw(20)<<std::real(coh);
            tfile_out<<std::setw(20)<<std::imag(coh);
            tfile_out<<std::setw(20)<<Jra[0];
            tfile_out<<std::setw(20)<<Jra[1];
            tfile_out<<std::setw(20)<<Jer[0];
            tfile_out<<std::setw(20)<<Jer[1];
            tfile_out<<std::endl;

            //prepare output streams
            /*std::string rho_t_fname=params.rhofile_fname;
            replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));
            std::ofstream rho_t_out(rho_t_fname);

            for(size_t ix=0; ix<params.Nx; ix++){
                double x=xgrid[ix];
                for(size_t iy=0; iy<params.Ny; iy++){
                    double y=ygrid[iy];



                    rho_t_out<<res(ix,iy)<<std::endl;
                }
            }

            rho_t_out.close();*/
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
