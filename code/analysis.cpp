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

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

        //prepare grids
        auto tgrid=create_grid(params.tmin,params.tmax,params.Nt);

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

        //create Dirac points
        std::vector<double> Dirac_Kx;
        std::vector<double> Dirac_Ky;
        std::vector<int> Dirac_type;

        gm.get_Dirac_points(Dirac_Kx,Dirac_Ky,Dirac_type);
        size_t nK=Dirac_Kx.size();

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

        for(size_t it=0; it<params.Nt; it++){
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
            }//loop over Dirac points
            std::cout<<std::endl;

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
