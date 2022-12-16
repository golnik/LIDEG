#include <iostream>
#include <vector>

#include "mini/ini.h"

#include "utils/utils.hpp"
#include "utils/multiarray.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "graphenemodel.hpp"
#include "model1/graphene.hpp"
#include "model1/WFs.hpp"
#include "model2/graphene2.hpp"
#include "Nlayer/nlayer.hpp"

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
        HexagonalTBModel* tb=new HexagonalTBModel(params.a);
        Graphene* gm;

        //create kgrid
        Grid2D* kxygrid;

        if(params.kgrid_type==kgrid_types::ucell){
            double Ox=0.;
            double Oy=0.;
            double b1x=2.*M_PI/(sqrt(3.)*params.a);
            double b1y=2.*M_PI/(params.a);
            double b2x=b1x;
            double b2y=-b1y;

            kxygrid=new UCellGrid2D(Ox,Oy,b1x,b1y,params.Nkx,b2x,b2y,params.Nky);
        }
        else{
            throw std::string("Integration is possible only for kgrid_type=ucell!");
        }

        //area of BZ
        Integrator2D* integrator=new Integrator2D(kxygrid);
        double SBZ=pow(2.*M_PI,2.)/(0.5*sqrt(3.)*params.a*params.a);

        if(params.model==models::hommelhoff){
            gm=new GrapheneModel(params.a,params.e2p,params.gamma,params.s,params.Td,Efield_x,Efield_y);
        }
        else if(params.model==models::nlayer){            
            double eps[]={params.e2p,params.e2p};
            double g[]={params.gamma,0.39/au2eV};
            double s[]={params.s,0.0};

            gm=new NGraphene(tb,params.nlayers,
                        kxygrid,
                        eps,g,s,
                        Efield_x,Efield_y);
        }

        size_t Nstates=gm->nstates();

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

            //integrate band populations
            integrator->trapz(
                [dens_vv](const size_t& ix, const size_t& iy){
                    return dens_vv(ix,iy);
                },
                pop0
            );

            integrator->trapz(
                [dens_cc](const size_t& ix, const size_t& iy){
                    return dens_cc(ix,iy);
                },
                pop1
            );                

            double coh_re=0.;
            double coh_im=0.;

            //integrate coherences
            integrator->trapz(
                [dens_cv_re](const size_t& ix, const size_t& iy){
                    return dens_cv_re(ix,iy);
                },
                coh_re
            );

            integrator->trapz(
                [dens_cv_im](const size_t& ix, const size_t& iy){
                    return dens_cv_im(ix,iy);
                },
                coh_im
            );

            coh+=coh_re+I*coh_im;

            /*for(int dir=0; dir<2; dir++){
                Jra[dir]+=kxygrid->integrate(
                    [dens_vv,dens_cc,
                    kxygrid,
                    time,
                    Afield_x,Afield_y,
                    gm,
                    dir](const size_t& ikx, const size_t& iky){
                        double kx0=(*kxygrid)(ikx,iky)[0];
                        double ky0=(*kxygrid)(ikx,iky)[1];

                        double kxt=kx0+(*Afield_x)(time);
                        double kyt=ky0+(*Afield_y)(time);

                        if(dir==0){
                            return 0.;
                            //return gm.px_vv(kxt,kyt)*dens_vv(ikx,iky)
                            //      +gm.px_cc(kxt,kyt)*dens_cc(ikx,iky);
                        }
                        else if(dir==1){
                            return 0.;
                            //return gm.py_vv(kxt,kyt)*dens_vv(ikx,iky)
                            //      +gm.py_cc(kxt,kyt)*dens_cc(ikx,iky);
                        }
                    }
                );
            }

            for(int dir=0; dir<2; dir++){
                Jer[dir]+=kxygrid->integrate(
                    [dens_cv_re,dens_cv_im,
                    kxygrid,
                    time,
                    Afield_x,Afield_y,
                    gm,
                    dir](const size_t& ikx, const size_t& iky){
                        double kx0=(*kxygrid)(ikx,iky)[0];
                        double ky0=(*kxygrid)(ikx,iky)[1];

                        double kxt=kx0+(*Afield_x)(time);
                        double kyt=ky0+(*Afield_y)(time);

                        complex_t dens_cv=dens_cv_re(ikx,iky)+I*dens_cv_im(ikx,iky);

                        if(dir==0){
                            return 0.;
                            //return 2.*std::real(gm.px_cv(kxt,kyt)*dens_cv);
                        }
                        else if(dir==1){
                            return 0.;
                            //return 2.*std::real(gm.py_cv(kxt,kyt)*dens_cv);
                        }
                    }
                );
            }*/

            tfile_out<<std::setw(20)<<time*au2fs;
            tfile_out<<std::setw(20)<<(*Afield_x)(time);
            tfile_out<<std::setw(20)<<(*Afield_y)(time);
            tfile_out<<std::setw(20)<<(*Efield_x)(time);
            tfile_out<<std::setw(20)<<(*Efield_y)(time);
            tfile_out<<std::setw(20)<<pop0/SBZ;
            tfile_out<<std::setw(20)<<pop1/SBZ;
            tfile_out<<std::setw(20)<<std::real(coh)/SBZ;
            tfile_out<<std::setw(20)<<std::imag(coh)/SBZ;
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
