#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>

#include "mini/ini.h"

#include "utils/utils.hpp"
#include "utils/grid.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "graphenemodel.hpp"
#include "model1/graphene.hpp"
//#include "model1/WFs.hpp"
#include "model2/graphene2.hpp"
#include "Nlayer/nlayer.hpp"

#include "uBLASBoostOdeint.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

#include <boost/format.hpp>

//#include <boost/filesystem.hpp>
//namespace fs=boost::filesystem;
#include <filesystem>
namespace fs=std::filesystem;

typedef matrix<complex_t> state_type;

typedef runge_kutta4<state_type> rk4;
typedef runge_kutta_cash_karp54<state_type> rkck54;
typedef runge_kutta_fehlberg78<state_type> rkf78;
typedef bulirsch_stoer<state_type> bst;
typedef runge_kutta_dopri5<state_type> rkd5;

typedef controlled_runge_kutta<rkck54> ctrl_rkck54;

typedef adams_bashforth_moulton<7,state_type> abm;

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

        //create output directory if does not exist
        fs::path outpath(params.outdir);
        fs::create_directory(outpath);

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

        HexagonalTBModel* tb=new HexagonalTBModel(params.a);
        Graphene* gm;

        //create kgrid
        Grid2D* kxygrid;

        if(params.kgrid_type==kgrid_types::kpoint){
            //create K and Kp points
            std::vector<double> Dirac_Kx;
            std::vector<double> Dirac_Ky;
            std::vector<int> Dirac_type;

            tb->get_Dirac_points(Dirac_Kx,Dirac_Ky,Dirac_type);

            //we simulate dynamics at K point only
            double Kx=Dirac_Kx[0];
            double Ky=Dirac_Ky[0];

            Grid1D* kx_grid=new RegularGrid1D(Kx-params.dkx,Kx+params.dkx,params.Nkx);
            Grid1D* ky_grid=new RegularGrid1D(Ky-params.dky,Ky+params.dky,params.Nky);

            kxygrid=new RegularGrid2D(kx_grid,ky_grid);
        }
        else if(params.kgrid_type==kgrid_types::ucell){
            double Ox=0.;
            double Oy=0.;
            double b1x=2.*M_PI/(sqrt(3.)*params.a);
            double b1y=2.*M_PI/(params.a);
            double b2x=b1x;
            double b2y=-b1y;

            kxygrid=new UCellGrid2D(Ox,Oy,b1x,b1y,params.Nkx,b2x,b2y,params.Nky);
        }

        if(params.model==models::hommelhoff){
            double e2p=params.e2p[0];
            double gamma=params.gamma[0];
            double s=params.s[0];
            gm=new GrapheneModel(params.a,e2p,gamma,s,params.Td,Efield_x,Efield_y);
        }
        else if(params.model==models::nlayer){
            gm=new NGraphene(tb,params.nlayers,
                        kxygrid,
                        params.e2p,params.gamma,params.s,params.Td,
                        Efield_x,Efield_y);
        }

        size_t Nstates=gm->nstates();

        //write grids to file
        std::ofstream grid_out(params.kgfile_fname);
        grid_out<<params.Nkx*params.Nky<<std::endl;

        grid_out<<std::scientific;
        grid_out<<std::setprecision(8);

        for(size_t ikx=0; ikx<params.Nkx; ikx++){
            for(size_t iky=0; iky<params.Nky; iky++){
                grid_out<<std::setw(20)<<(*kxygrid)(ikx,iky)[0];
                grid_out<<std::setw(20)<<(*kxygrid)(ikx,iky)[1]<<std::endl;
            }
        }
        grid_out.close();

        //prepare initial densities
        if (params.T <= 0)
        {
            matrix<state_type> rho_t_kxky(params.Nkx, params.Nky);
            for (size_t ikx = 0; ikx < params.Nkx; ikx++)
            {
                for (size_t iky = 0; iky < params.Nky; iky++)
                {
                    rho_t_kxky(ikx, iky) = state_type(Nstates, Nstates);
                    for (size_t i = 0; i < Nstates; i++)
                    {
                        for (size_t j = 0; j < Nstates; j++)
                        {
                            rho_t_kxky(ikx, iky)(i, j) = 0.;
                        }
                    }
                    rho_t_kxky(ikx, iky)(0, 0) = 1.;
                }
            }
        }
        else
        {
            matrix<state_type> rho_t_kxky(params.Nkx, params.Nky);
            for (size_t ikx = 0; ikx < params.Nkx; ikx++)
            {
                for (size_t iky = 0; iky < params.Nky; iky++)
                {
                    double kx0 = (*kxygrid)(ikx, iky)[0];
                    double ky0 = (*kxygrid)(ikx, iky)[1];
                    double kBT = KB * params.T;
                    rho_t_kxky(ikx, iky) = state_type(Nstates, Nstates);
                    for (size_t i = 0; i < Nstates; i++)
                    {
                        for (size_t j = 0; j < Nstates; j++)
                        {
                            if (i == j && i == 1)
                            {
                                double fc = 1 / (exp(gm->get_energy(kx0, ky0, i) / kBT) + 1);
                                rho_t_kxky(ikx, iky)(i, j) = fc;
                            }
                            else
                            {
                                rho_t_kxky(ikx, iky)(i, j) = 0.;
                            }
                        }
                    }
                    double fv = 1 / (exp(gm->get_energy(kx0, ky0, 0) / kBT) + 1);
                    rho_t_kxky(ikx, iky)(0, 0) = fv;
                }
            }
        }

        auto tgrid=create_grid(params.tmin,params.tmax,params.Nt);
        double dt=(params.tmax-params.tmin)/(params.Nt-1);
        for(size_t it=0; it<params.Nt; it++){//time loop
            //output rho for every timestep
            std::string dens_t_fname=params.densfile_fname;
            replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));
            std::ofstream dens_t_out(dens_t_fname);

            std::cout<<"Time step "<<it+1<<" out of "<<params.Nt<<std::endl;

            double time=tgrid[it];

            #pragma omp parallel for
            for(size_t ikx=0; ikx<params.Nkx; ikx++){
                for(size_t iky=0; iky<params.Nky; iky++){
                    //double kx0=(*kx_grid)[ikx];
                    //double ky0=(*ky_grid)[iky];

                    double kx0=(*kxygrid)(ikx,iky)[0];
                    double ky0=(*kxygrid)(ikx,iky)[1];

                    auto system=[&gm,kx0,ky0,Afield_x,Afield_y](const state_type& rho, state_type& drhodt, const double t){
                        double kxt=kx0+(*Afield_x)(t);
                        double kyt=ky0+(*Afield_y)(t);

                        gm->propagate(rho,drhodt,t,kxt,kyt);
                        //gm2.propagate(rho,drhodt,t,kx0,ky0);
                        //ngm.propagate(rho,drhodt,t,kxt,kyt);
                        return;
                    };

                    integrate_adaptive(
                    //integrate_const(
                        //rk4(),
                        //rkck54(),
                        //ctrl_rkck54(),
                        //rkf78(),
                        //rkd5(),
                        //abm(),
                        make_controlled(params.err_abs,params.err_rel,rkf78()),
                        system,rho_t_kxky(ikx,iky),time,time+dt,params.ddt);
                }
            }

            //write kgrid files
            dens_t_out<<std::scientific;
            //rho_t_out<<std::fixed;
            dens_t_out<<std::setprecision(12);

            //write file header
            dens_t_out<<"#";
            
            size_t col=1;
            for(size_t ist=0; ist<Nstates; ist++){
                std::string pop_str="dens["+std::to_string(ist+1)+"]("+std::to_string(col)+")";
                dens_t_out<<std::setw(20)<<pop_str;
                col++;
            }

            for(size_t ist=0; ist<Nstates; ist++){
                for(size_t jst=ist+1; jst<Nstates; jst++){
                    std::string coh_re_str="Re{coh["+std::to_string(ist+1)+","+std::to_string(jst+1)+"]("+std::to_string(col)+")}";
                    col++;
                    std::string coh_im_str="Im{coh["+std::to_string(ist+1)+","+std::to_string(jst+1)+"]("+std::to_string(col)+")}";
                    col++;

                    dens_t_out<<std::setw(20)<<coh_re_str;
                    dens_t_out<<std::setw(20)<<coh_im_str;
                }
            }

            //dfdk for conductivity calculations
            for (size_t ist = 0; ist < Nstates; ist++)
            {
                std::string pop_str = "dfdk[" + std::to_string(ist + 1) + "](" + std::to_string(col) + ")";
                dens_t_out << std::setw(20) << pop_str;
                col++;
            }

            for (size_t ist = 0; ist < Nstates; ist++)
            {
                for (size_t jst = ist + 1; jst < Nstates; jst++)
                {
                    std::string coh_re_str = "Re{dfdk[" + std::to_string(ist + 1) + "," + std::to_string(jst + 1) + "](" + std::to_string(col) + ")}";
                    col++;
                    std::string coh_im_str = "Im{dfdk[" + std::to_string(ist + 1) + "," + std::to_string(jst + 1) + "](" + std::to_string(col) + ")}";
                    col++;

                    dens_t_out << std::setw(20) << coh_re_str;
                    dens_t_out << std::setw(20) << coh_im_str;
                }
            }
            dens_t_out<<std::endl;

            //write data
            for(size_t ikx=0; ikx<params.Nkx; ikx++){
                for(size_t iky=0; iky<params.Nky; iky++){
                    for(size_t ist=0; ist<Nstates; ist++){
                        double rho=std::abs(rho_t_kxky(ikx,iky)(ist,ist));
                        dens_t_out<<std::setw(20)<<rho;
                    }
                    for(size_t ist=0; ist<Nstates; ist++){
                        for(size_t jst=ist+1; jst<Nstates; jst++){
                            double re=std::real(rho_t_kxky(ikx,iky)(ist,jst));
                            double im=std::imag(rho_t_kxky(ikx,iky)(ist,jst));
                            dens_t_out<<std::setw(20)<<re;
                            dens_t_out<<std::setw(20)<<im;
                        }
                    }

                    /////////////////

                    /// Cal df/dk

                    /////////////////

/*                     double kx0=(*kxygrid)(ikx,iky)[0];
                    double ky0=(*kxygrid)(ikx,iky)[1];

                    double kxt=kx0+(*Afield_x)(time);
                    double kyt=ky0+(*Afield_y)(time);

                    auto const kx_ad=make_fvar<double,1>(kxt);
                    auto const autodiff=std::abs(rho_t_kxky(ikx, iky)(1, 1));
                    double df_dkx = autodiff.derivative(1); */

                    //df_dk output to data file
                    //forward difference formula (modify to central difference to check the accuracy)
                    if (ikx < params.Nkx - 1)
                    {
                        double df = std::abs(rho_t_kxky(ikx + 1, iky)(1, 1)) - std::abs(rho_t_kxky(ikx, iky)(1, 1));
                        // double dk = (*kxygrid)(ikx + 1, iky)[0] - (*kxygrid)(ikx, iky)[0];
                        double dk = 2. * M_PI / (sqrt(3.) * params.a * params.Nkx);
                        dens_t_out << std::setw(20) << -df / dk;
                        dens_t_out << std::setw(20) << df / dk;
                    }
                    else
                    {
                        dens_t_out << std::setw(20) << 1e-15;
                        dens_t_out << std::setw(20) << 1e-15;
                    }

                    //off-diagonal elements of df/dk
                    for (size_t ist = 0; ist < Nstates; ist++)
                    {
                        for (size_t jst = ist + 1; jst < Nstates; jst++)
                        {
                            if (ikx < params.Nkx - 1)
                            {
                                double dk = 2. * M_PI / (sqrt(3.) * params.a * params.Nkx);
                                double dfre = std::real(rho_t_kxky(ikx + 1, iky)(ist, jst)) - std::real(rho_t_kxky(ikx, iky)(ist, jst));
                                double dfim = std::imag(rho_t_kxky(ikx + 1, iky)(ist, jst)) - std::imag(rho_t_kxky(ikx, iky)(ist, jst));
                                dens_t_out << std::setw(20) << dfre / dk;
                                dens_t_out << std::setw(20) << dfim / dk;
                            }
                            else
                            {
                                dens_t_out << std::setw(20) << 1e-15;
                                dens_t_out << std::setw(20) << 1e-15;
                            }
                        }
                    }

                    dens_t_out<<std::endl;
                }
            }

            //close streams
            dens_t_out.close();
        }//time loop
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
