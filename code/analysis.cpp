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

        //prepare grids
        auto tgrid_tmp=create_grid(params.tmin,params.tmax,params.Nt);
        std::vector<double> tgrid;

        //std::vector<int> tsteps(params.Nt);
        //std::iota(tsteps.begin(),tsteps.end(),0);

        std::vector<int> tsteps={tstep};
        
        for(auto it: tsteps){
            tgrid.push_back(tgrid_tmp[it]);
        }

        //std::vector<double> tgrid=

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

        //create graphene layer
        double R0x=0.5*(params.xmax+params.xmin);
        double R0y=0.5*(params.ymax+params.ymin);        
        GrapheneLayer gl(params.a,
                params.Nclx,params.Ncly,
                R0x,R0y,params.Rmax);

        //area of BZ
        double SBZ=pow(2.*M_PI,2.)/(0.5*sqrt(3.)*params.a*params.a);

        //create Dirac points
        std::vector<double> Dirac_Kx;
        std::vector<double> Dirac_Ky;
        std::vector<int> Dirac_type;

        gm.get_Dirac_points(Dirac_Kx,Dirac_Ky,Dirac_type);
        size_t nK=Dirac_Kx.size();

        std::vector<WFs_grid*> wfs(nK);

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

            //WFs wfs(&gm,&gl,params.Z);
            //WFs_grid wfs_g(&gm,&gl,params.Z,kx_grid,ky_grid);
            //wfs.push_back(WFs_grid(&gm,&gl,params.Z,kx_grid,ky_grid));
            wfs[iK]=new WFs_grid(&gm,&gl,params.Z,kx_grid,ky_grid);
        }

        //print atom positions to file
        std::ofstream atoms_out;
        atoms_out.open(params.afile_fname);
        gl.print_atoms(atoms_out);
        atoms_out.close();

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
        tfile_out<<std::setw(20)<<"dens_vv";
        tfile_out<<std::setw(20)<<"dens_cc";
        tfile_out<<std::setw(20)<<"Re{dens_cv}";
        tfile_out<<std::setw(20)<<"Im{dens_cv}";
        tfile_out<<std::setw(20)<<"Jx_intra";
        tfile_out<<std::setw(20)<<"Jy_intra";
        tfile_out<<std::setw(20)<<"Jx_inter";
        tfile_out<<std::setw(20)<<"Jy_inter";
        tfile_out<<std::endl;

        /*MultiArray<complex_t,64,64,64,64,64> tmp;

        for(size_t ikx=0; ikx<params.Nkx; ikx++){
            for(size_t iky=0; iky<params.Nky; iky++){
                for(size_t ix=0; ix<params.Nx; ix++){
                    for(size_t iy=0; iy<params.Ny; iy++){
                        for(size_t iz=0; iz<params.Nz; iz++){
                            double kx=
                        }
                    }
                }
            }
        }

        return 1;*/

        //array for real space data
        //matrix<double> res(params.Nx,params.Ny);
        MultiArray<double,64,64,64> res;

        for(auto it: tsteps){
            //read density from file
            std::string dens_t_fname=params.densfile_fname;
            replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));

            std::cout<<dens_t_fname<<std::endl;

            matrix_t dens_vv=read_2D_from_file<matrix_t>(dens_t_fname,0,params.Nkx,params.Nky);
            matrix_t dens_cc=read_2D_from_file<matrix_t>(dens_t_fname,1,params.Nkx,params.Nky);
            matrix_t dens_cv_re=read_2D_from_file<matrix_t>(dens_t_fname,2,params.Nkx,params.Nky);
            matrix_t dens_cv_im=read_2D_from_file<matrix_t>(dens_t_fname,3,params.Nkx,params.Nky);

            double time=tgrid[it];

            //integrate over kgrids
            double pop0=0.;
            double pop1=0.;
            complex_t coh=0.;
            double Jra[2]={0.,0.};
            double Jer[2]={0.,0.};

            double* mres=new double [Nmulti];

            for(size_t ix=0; ix<params.Nx; ix++){
                for(size_t iy=0; iy<params.Ny; iy++){
                    for(size_t iz=0; iz<params.Nz; iz++){
                        res(ix,iy,iz)=0.;
                    }
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

                /*double kxmin=-2.;
                double kxmax=2.;
                double kymin=-2.;
                double kymax=2.;*/

                auto kx_grid=create_grid(kxmin,kxmax,params.Nkx,params.kgrid_type);
                auto ky_grid=create_grid(kymin,kymax,params.Nky,params.kgrid_type);

                //create WFs model
                //WFs wfs(&gm,&gl,params.Z);
                //WFs_grid wfs_g(&gm,&gl,params.Z,kx_grid,ky_grid);

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
                    [dens_vv](const size_t& ix, const size_t& iy){
                        return dens_vv(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax
                );

                pop1+=integrate(params.Nkx,params.Nky,
                    [dens_cc](const size_t& ix, const size_t& iy){
                        return dens_cc(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax                    
                );                

                //integrate coherences
                double coh_re=integrate(params.Nkx,params.Nky,
                    [dens_cv_re](const size_t& ix, const size_t& iy){
                        return dens_cv_re(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax
                );

                double coh_im=integrate(params.Nkx,params.Nky,
                    [dens_cv_im](const size_t& ix, const size_t& iy){
                        return dens_cv_im(ix,iy);
                    },
                    kxmin,kxmax,
                    kymin,kymax
                );

                coh+=coh_re+I*coh_im;

                for(int dir=0; dir<2; dir++){
                    Jra[dir]+=integrate(params.Nkx,params.Nky,
                        [dens_vv,dens_cc,
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
                                return gm.px_vv(kxt,kyt)*dens_vv(ikx,iky)
                                      +gm.px_cc(kxt,kyt)*dens_cc(ikx,iky);
                            }
                            else if(dir==1){
                                return gm.py_vv(kxt,kyt)*dens_vv(ikx,iky)
                                      +gm.py_cc(kxt,kyt)*dens_cc(ikx,iky);
                            }
                        },
                        kxmin,kxmax,
                        kymin,kymax
                    );
                }

                for(int dir=0; dir<2; dir++){
                    Jer[dir]+=integrate(params.Nkx,params.Nky,
                        [dens_cv_re,dens_cv_im,
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

                            complex_t dens_cv=dens_cv_re(ikx,iky)+I*dens_cv_im(ikx,iky);

                            if(dir==0){
                                return 2.*std::real(gm.px_cv(kxt,kyt)*dens_cv);
                            }
                            else if(dir==1){
                                return 2.*std::real(gm.py_cv(kxt,kyt)*dens_cv);
                            }
                        },
                        kxmin,kxmax,
                        kymin,kymax
                    );
                }

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
                                size_t ikx=mgrid[im];
                                size_t iky=mgrid[im+Nmulti];

                                double kx=kx_grid[ikx];
                                double ky=ky_grid[iky];

                                double kxt=kx+(*Afield_x)(time);
                                double kyt=ky+(*Afield_y)(time);

                                complex_t dens_cv=dens_cv_re(ikx,iky)+I*dens_cv_im(ikx,iky);

                                //mres[im]=wfs.phi_2pz(x,y,z);

                                //mres[im]=wfs.PhiA1(x,y,z,kx,ky);
                                //mres[im]=wfs.PhiA2(x,y,z,kx,ky);

                                //double rho_vv=std::norm(wfs.psip(x,y,z,kx,ky));
                                //double rho_cc=std::norm(wfs.psim(x,y,z,kx,ky));

                                //complex_t psip=wfs.psip(x,y,z,kxt,kyt);
                                //complex_t psim=wfs.psim(x,y,z,kxt,kyt);

                                //complex_t rho_vc=std::conj(wfs.psip(x,y,z,kx,ky))*wfs.psim(x,y,z,kx,ky);

                                complex_t psip=wfs[iK]->psip(x,y,z,ikx,iky);
                                complex_t psim=wfs[iK]->psim(x,y,z,ikx,iky);

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
            std::cout<<std::endl;

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

            //output real space data
            std::string rho_t_fname=params.rhofile_fname;
            replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));
            std::ofstream rho_t_out(rho_t_fname);

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
