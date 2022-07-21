#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>

#include "mini/ini.h"

#include "utils.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "model1/graphene.hpp"
#include "model1/WFs.hpp"

#include "uBLASBoostOdeint.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

#include <boost/format.hpp>

typedef matrix<complex_t> state_type;

typedef runge_kutta4<state_type> rk4;
typedef runge_kutta_cash_karp54<state_type> rkck54;
typedef runge_kutta_fehlberg78<state_type> rkf78;
typedef bulirsch_stoer<state_type> bst;
typedef runge_kutta_dopri5<state_type> rkd5;

typedef controlled_runge_kutta<rkck54> ctrl_rkck54;

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

        //prepare grids
        auto tgrid  =create_grid(params.tmin,params.tmax,params.Nt);
        auto kx_grid=create_grid(params.kx_min,params.kx_max,params.Nkx);
        auto ky_grid=create_grid(params.ky_min,params.ky_max,params.Nky);

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

        //prepare initial densities
        matrix<state_type> rho_t_kxky(params.Nkx,params.Nky);
        for(size_t ikx=0; ikx<params.Nkx; ikx++){
            for(size_t iky=0; iky<params.Nky; iky++){
                rho_t_kxky(ikx,iky)=state_type(2,2);
                for(size_t i=0; i<2; i++){
                    for(size_t j=0; j<2; j++){
                        rho_t_kxky(ikx,iky)(i,j)=0.;
                    }
                }
                rho_t_kxky(ikx,iky)(0,0)=1.;
            }
        }

        //prepare output streams
        std::ofstream tfile_out;
        tfile_out.open(params.tfile_fname);

        double dt=(params.tmax-params.tmin)/(params.Nt-1);
        for(size_t it=0; it<params.Nt; it++){
            //output rho for every timestep
            std::string rho_t_fname=params.rhofile_fname;
            replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));
            std::ofstream rho_t_out(rho_t_fname);

            std::cout<<"Time step "<<it+1<<" out of "<<params.Nt<<std::endl;

            double time=tgrid[it];

            #pragma omp parallel for
            for(size_t ikx=0; ikx<params.Nkx; ikx++){
                for(size_t iky=0; iky<params.Nky; iky++){
                    double kx0=kx_grid[ikx];
                    double ky0=ky_grid[iky];

                    auto system=[gm,kx0,ky0,Afield_x,Afield_y](const state_type& rho, state_type& drhodt, const double t){
                        double kxt=kx0+(*Afield_x)(t);
                        double kyt=ky0+(*Afield_y)(t);

                        gm.propagate(rho,drhodt,t,kxt,kyt);
                        return;
                    };

                    integrate_adaptive(
                    //integrate_const(
                        //rk4(),
                        //rkck54(),
                        //ctrl_rkck54(),
                        //rkf78(),
                        //rkd5(),
                        make_controlled(params.err_abs,params.err_rel,rkd5()), 
                        system,rho_t_kxky(ikx,iky),time,time+dt,params.ddt);
                }
            }

            double pop0=0.;
            double pop1=0.;
            complex_t coh=0.;

            for(size_t ikx=0; ikx<params.Nkx; ikx++){
                for(size_t iky=0; iky<params.Nky; iky++){
                    pop0+=abs(rho_t_kxky(ikx,iky)(0,0));
                    pop1+=abs(rho_t_kxky(ikx,iky)(1,1));
                    coh+=rho_t_kxky(ikx,iky)(0,1);

                    //write kgrid files
                    //rho_t_out<<kx_grid[ikx]<<" "<<ky_grid[iky]<<" ";
                    rho_t_out<<std::abs(rho_t_kxky(ikx,iky)(0,0))<<" ";
                    rho_t_out<<std::abs(rho_t_kxky(ikx,iky)(1,1))<<" ";
                    rho_t_out<<std::real(rho_t_kxky(ikx,iky)(0,1))<<" ";
                    rho_t_out<<std::imag(rho_t_kxky(ikx,iky)(0,1))<<" ";
                    rho_t_out<<std::endl;
                }
            }

            tfile_out<<std::fixed;
            tfile_out<<std::setprecision(5);
            tfile_out<<time*au2fs<<" ";
            tfile_out<<(*Afield_x)(time)<<" ";
            tfile_out<<(*Afield_y)(time)<<" ";
            tfile_out<<(*Efield_x)(time)<<" ";
            tfile_out<<(*Efield_y)(time)<<" ";
            tfile_out<<pop0<<" ";
            tfile_out<<pop1<<" ";
            tfile_out<<std::real(coh)<<" ";
            tfile_out<<std::imag(coh)<<" ";

            tfile_out<<std::endl;

            //close streams
            rho_t_out.close();
        }

        //close streams
        tfile_out.close();
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
