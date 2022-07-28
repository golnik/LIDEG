#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>

#include "mini/ini.h"

#include "utils/utils.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "model1/graphene.hpp"
#include "model1/WFs.hpp"
#include "model2/graphene2.hpp"

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

        GrapheneModel2 gm2(params.a,params.e2p,params.gamma,params.s,
                           Efield_x,Efield_y,
                           Afield_x,Afield_y);

        //create K and Kp points
        std::vector<double> Dirac_Kx;
        std::vector<double> Dirac_Ky;
        std::vector<int> Dirac_type;

        gm.get_Dirac_points(Dirac_Kx,Dirac_Ky,Dirac_type);

        //we simulate dynamics at K point only
        double Kx=Dirac_Kx[0];
        double Ky=Dirac_Ky[0];

        //create k grids
        auto kx_grid=create_grid(Kx-params.dkx,Kx+params.dkx,params.Nkx,params.kgrid_type);
        auto ky_grid=create_grid(Ky-params.dky,Ky+params.dky,params.Nky,params.kgrid_type);

        //write grids to file
        std::ofstream grid_out("output/grids.dat");
        //write kx grid
        grid_out<<kx_grid.size()<<std::endl;
        for(size_t ikx=0; ikx<params.Nkx; ikx++)
            grid_out<<kx_grid[ikx]-Kx<<" ";
        grid_out<<std::endl;
        grid_out<<ky_grid.size()<<std::endl;
        for(size_t iky=0; iky<params.Nky; iky++)
            grid_out<<ky_grid[iky]-Ky<<" ";
        grid_out<<std::endl;  
        grid_out.close();

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
        
        std::ofstream test_out("test.out");

        auto tgrid=create_grid(params.tmin,params.tmax,params.Nt);
        double dt=(params.tmax-params.tmin)/(params.Nt-1);
        for(size_t it=0; it<params.Nt; it++){//time loop
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

                    auto system=[gm,gm2,kx0,ky0,Afield_x,Afield_y](const state_type& rho, state_type& drhodt, const double t){
                        double kxt=kx0+(*Afield_x)(t);
                        double kyt=ky0+(*Afield_y)(t);

                        gm.propagate(rho,drhodt,t,kxt,kyt);
                        //gm2.propagate(rho,drhodt,t,kx0,ky0);
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

            //write kgrid files
            rho_t_out<<std::scientific;
            //rho_t_out<<std::fixed;
            rho_t_out<<std::setprecision(8);
            rho_t_out<<"#"<<std::setw(19)<<"rho_vv";
            rho_t_out<<     std::setw(20)<<"rho_cc";
            rho_t_out<<     std::setw(20)<<"Re{rho_cv}";
            rho_t_out<<     std::setw(20)<<"Im{rho_cv}";
            rho_t_out<<std::endl;
            for(size_t ikx=0; ikx<params.Nkx; ikx++){
                for(size_t iky=0; iky<params.Nky; iky++){
                    double rho_vv=std::abs(rho_t_kxky(ikx,iky)(0,0));
                    double rho_cc=std::abs(rho_t_kxky(ikx,iky)(1,1));
                    complex_t rho_cv=rho_t_kxky(ikx,iky)(0,1);

                    /*double rho_vv=std::norm(rho_t_kxky(ikx,iky)(0,0));
                    double rho_cc=std::norm(rho_t_kxky(ikx,iky)(1,1));
                    complex_t rho_cv=std::conj(rho_t_kxky(ikx,iky)(0,0))
                                    *rho_t_kxky(ikx,iky)(1,1);*/

                    //rho_t_out<<kx_grid[ikx]<<" "<<ky_grid[iky]<<" ";
                    rho_t_out<<std::setw(20)<<rho_vv;
                    rho_t_out<<std::setw(20)<<rho_cc;
                    rho_t_out<<std::setw(20)<<std::real(rho_cv);
                    rho_t_out<<std::setw(20)<<std::imag(rho_cv);
                    rho_t_out<<std::setw(20)<<std::abs(rho_cv);
                    rho_t_out<<std::setw(20)<<std::arg(rho_cv);

                    rho_t_out<<std::endl;
                }
            }

            test_out<<std::scientific;
            test_out<<std::setprecision(8);
            test_out<<time<<" "<<std::abs(rho_t_kxky(0,0)(0,0))<<std::endl;

            //close streams
            rho_t_out.close();
        }//time loop

        test_out.close();
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
