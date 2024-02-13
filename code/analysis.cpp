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

        //prepare output streams
        std::ofstream tfile_out;
        tfile_out.open(params.tfile_fname);

        tfile_out<<std::scientific;
        //tfile_out<<std::fixed;
        tfile_out<<std::setprecision(8);

        tfile_out<<"#";
        tfile_out<<std::setw(19)<<"Time(1)";
        tfile_out<<std::setw(20)<<"Ax(2)";
        tfile_out<<std::setw(20)<<"Ay(3)";
        tfile_out<<std::setw(20)<<"Ex(4)";
        tfile_out<<std::setw(20)<<"Ey(5)";

        tfile_out<<std::setw(20)<<"Jx_intra(6)";
        tfile_out<<std::setw(20)<<"Jy_intra(7)";
        tfile_out<<std::setw(20)<<"Jx_inter(8)";
        tfile_out<<std::setw(20)<<"Jy_inter(9)";

        size_t col=10;
        for(size_t ist=0; ist<Nstates; ist++){
            std::string pop_str="dens["+std::to_string(ist+1)+"]("+std::to_string(col)+")";
            tfile_out<<std::setw(20)<<pop_str;
            col++;
        }

        for(size_t ist=0; ist<Nstates; ist++){
            for(size_t jst=ist+1; jst<Nstates; jst++){
                std::string coh_re_str="coh_re["+std::to_string(ist+1)+","+std::to_string(jst+1)+"]("+std::to_string(col)+")";
                col++;
                std::string coh_im_str="coh_im["+std::to_string(ist+1)+","+std::to_string(jst+1)+"]("+std::to_string(col)+")";
                col++;

                tfile_out<<std::setw(20)<<coh_re_str;
                tfile_out<<std::setw(20)<<coh_im_str;
            }
        }

        // conductivity (numbering of columns should be modified)
        tfile_out << std::setw(20) << "Re{sig}(14)";
        tfile_out << std::setw(20) << "Im{sig}(15)";
        tfile_out << std::setw(20) << "Abs{sig}(16)";

        tfile_out<<std::endl;

        std::vector<matrix_t> dens_data(Nstates);
        std::vector<matrix_t> coh_re_data(Nstates*(Nstates-1)/2);
        std::vector<matrix_t> coh_im_data(Nstates*(Nstates-1)/2);
        // store df/dk
        std::vector<matrix_t> df_dk(Nstates);

        std::vector<double> pops(Nstates);
        std::vector<complex_t> cohs(Nstates*(Nstates-1)/2);

        for(size_t it=0; it<params.Nt; it++){
            double time=tgrid[it];

            //read density from file
            std::string dens_t_fname=params.densfile_fname;
            replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));

            std::cout<<dens_t_fname<<std::endl;

            size_t col=0;
            for(size_t ist=0; ist<Nstates; ist++){
                dens_data[ist]=read_2D_from_file<matrix_t>(dens_t_fname,col,params.Nkx,params.Nky);
                col++;
            }

            size_t indx=0;
            for(size_t ist=0; ist<Nstates; ist++){
                for(size_t jst=ist+1; jst<Nstates; jst++){
                    coh_re_data[indx]=read_2D_from_file<matrix_t>(dens_t_fname,col,params.Nkx,params.Nky);
                    col++;
                    coh_im_data[indx]=read_2D_from_file<matrix_t>(dens_t_fname,col,params.Nkx,params.Nky);
                    col++;
                    indx++;
                }
            }

            // read and save df/dk
            for (size_t ist = 0; ist < Nstates; ist++)
            {
                df_dk[ist] = read_2D_from_file<matrix_t>(dens_t_fname, col, params.Nkx, params.Nky);
                col++;
            }
            double Nk = params.Nkx;

            /////////////////////////

            /// Calculations

            /////////////////////////

            //integrate band populations
            for(size_t ist=0; ist<Nstates; ist++){
                pops[ist]=0.;
                integrator->trapz(
                    [ist,dens_data](const size_t& ikx, const size_t& iky){
                        return dens_data[ist](ikx,iky);
                    },
                    pops[ist]
                );
            }

            //integrate coherences
            indx=0;
            for(size_t ist=0; ist<Nstates; ist++){
                for(size_t jst=ist+1; jst<Nstates; jst++){
                    double coh_re=0.;
                    double coh_im=0.;

                    integrator->trapz(
                        [indx,coh_re_data](const size_t& ikx, const size_t& iky){
                            return coh_re_data[indx](ikx,iky);
                        },
                        coh_re
                    );

                    integrator->trapz(
                        [indx,coh_im_data](const size_t& ikx, const size_t& iky){
                            return coh_im_data[indx](ikx,iky);
                        },
                        coh_im
                    );

                    cohs[indx]=coh_re+I*coh_im;

                    indx++;
                }
            }

            //integrate intraband current
            double Jra[2]={0.,0.};
            for(int dir=0; dir<2; dir++){
                integrator->trapz(
                    [dens_data,
                    kxygrid,
                    time,
                    Afield_x,Afield_y,
                    gm,
                    dir,Nstates](const size_t& ikx, const size_t& iky){
                        double kx0=(*kxygrid)(ikx,iky)[0];
                        double ky0=(*kxygrid)(ikx,iky)[1];

                        double kxt=kx0+(*Afield_x)(time);
                        double kyt=ky0+(*Afield_y)(time);

                        double res=0.;
                        for(size_t ist=0; ist<Nstates; ist++){
                            res+=gm->get_energy_grad(kxt,kyt,ist,dir)*dens_data[ist](ikx,iky);
                        }

                        return res;
                    },
                    Jra[dir]
                );
            }

            //integrate interband current
            double Jer[2]={0.,0.};
            for(int dir=0; dir<2; dir++){
                integrator->trapz(
                    [coh_re_data,coh_im_data,
                    kxygrid,
                    time,
                    Afield_x,Afield_y,
                    gm,
                    dir,Nstates](const size_t& ikx, const size_t& iky){
                        double kx0=(*kxygrid)(ikx,iky)[0];
                        double ky0=(*kxygrid)(ikx,iky)[1];

                        double kxt=kx0+(*Afield_x)(time);
                        double kyt=ky0+(*Afield_y)(time);

                        double res=0.;

                        size_t indx=0;
                        for(size_t ist=0; ist<Nstates; ist++){
                            for(size_t jst=ist+1; jst<Nstates; jst++){
                                double re=coh_re_data[indx](ikx,iky);
                                double im=coh_im_data[indx](ikx,iky);
                                complex_t dens=re+I*im;

                                complex_t P=I*(gm->get_energy(kxt,kyt,ist)-gm->get_energy(kxt,kyt,jst))
                                    *gm->get_dipole(kxt,kyt,ist,jst,dir);

                                res+=2.*std::real(P*dens);

                                indx++;
                            }
                        }

                        return res;
                    },
                    Jer[dir]
                );
            }

            ////////////////

            /// Conducutivity

            ////////////////

            std::complex<double> sigma[4] = {0., 0., 0., 0.};
            for (int dir = 0; dir < 1; dir++)
            {
                integrator->trapz(
                    [coh_re_data, coh_im_data, dens_data, df_dk, Nk,
                     kxygrid,
                     time,
                     Afield_x, Afield_y,
                     gm,
                     dir, Nstates](const size_t &ikx, const size_t &iky)
                    {
                        double kx0 = (*kxygrid)(ikx, iky)[0];
                        double ky0 = (*kxygrid)(ikx, iky)[1];
                        double kxt = kx0 + (*Afield_x)(time);
                        double kyt = ky0 + (*Afield_y)(time);

/*                         std::complex<double> rho10 = coh_re_data[0](ikx, iky) - I * coh_im_data[0](ikx, iky);
                        std::complex<double> rho01 = coh_re_data[0](ikx, iky) + I * coh_im_data[0](ikx, iky);
                        double rho00 = dens_data[0](ikx, iky);
                        double rho11 = dens_data[1](ikx, iky);

                        std::complex<double> rho[2][2];
                        rho[0][0] = dens_data[0](ikx, iky);
                        rho[1][1] = dens_data[1](ikx, iky);
                        rho[1][0] = coh_re_data[0](ikx, iky) + I * coh_im_data[0](ikx, iky);
                        rho[0][1] = coh_re_data[0](ikx, iky) - I * coh_im_data[0](ikx, iky); */

                        double omega = 0;
                        double Eta = 1e-2;
                        std::complex<double> diag_res = 0.;
/*                         std::complex<double> total_res_1 = 0.;
                        std::complex<double> total_res_2 = 0.;
                        std::complex<double> result = 0.; */

                        std::complex<double> J[2][2];
                        std::complex<double> J_corr[2][2];

                        for (size_t m = 0; m < Nstates; m++)
                        {
                            for (size_t n = 0; n < Nstates; n++)
                            {
                                if (m == n)
                                {
                                    J[m][n] = gm->get_energy_grad(kxt, kyt, m, dir);
                                }
                                else if (m != n)
                                {
                                    J[m][n] = (gm->get_energy(kxt, kyt, m) - gm->get_energy(kxt, kyt, n)) * gm->get_dipole(kxt, kyt, m, n, 0);
                                }
                            }
                        }

                        //// Diag terms

                        for (size_t m = 0; m < Nstates; m++)
                        {
                            for (size_t n = 0; n < Nstates; n++)
                            {
                                if (m == n)
                                {
                                    // double kBT = 0.02659;
                                    // double anltic_dfde = -(exp(gm->get_energy(kxt, kyt, m)) / kBT) /
                                    //                      ((pow(1 + exp(gm->get_energy(kxt, kyt, m)) / kBT, 2.0)) * kBT);

                                    double df_de = df_dk[m](ikx, iky) / gm->get_energy_grad(kxt, kyt, m, dir);

/*                                     if (ikx == Nk && m == 1)
                                    {
                                        df_de = df_dk[m](Nk, iky) / gm->get_energy_grad(kxt, kyt, m, dir);
                                    }
                                    else if (ikx == Nk && m == 0)
                                    {
                                        df_de = df_dk[m](Nk, iky) / gm->get_energy_grad(kxt, kyt, m, dir);                                        
                                    } */

                                    if (gm->get_energy_grad(kxt, kyt, m, dir) == 0)
                                    {
                                        df_de = 0;
                                    }

                                    // df_de = 0; //for inter band conductivity

                                    diag_res += df_de * J[n][m] * J[m][n] /
                                                ((gm->get_energy(kxt, kyt, n) - gm->get_energy(kxt, kyt, m) + omega + I * Eta));
                                }
                                else if (m != n)
                                {
                                    diag_res += (dens_data[m](ikx, iky) - dens_data[n](ikx, iky)) * J[n][m] * J[m][n] /
                                                ((gm->get_energy(kxt, kyt, n) - gm->get_energy(kxt, kyt, m)) * (gm->get_energy(kxt, kyt, n) - gm->get_energy(kxt, kyt, m) + omega + I * Eta));
                                }
                            }
                        }

                        return diag_res / I;
                    },
                    sigma[0]);
            }

            ////////////////

            /// write data

            ////////////////


            tfile_out<<std::setw(20)<<time*au2fs;
            tfile_out<<std::setw(20)<<(*Afield_x)(time);
            tfile_out<<std::setw(20)<<(*Afield_y)(time);
            tfile_out<<std::setw(20)<<(*Efield_x)(time);
            tfile_out<<std::setw(20)<<(*Efield_y)(time);

            tfile_out<<std::setw(20)<<Jra[0];
            tfile_out<<std::setw(20)<<Jra[1];
            tfile_out<<std::setw(20)<<Jer[0];
            tfile_out<<std::setw(20)<<Jer[1];

            for(size_t ist=0; ist<Nstates; ist++){
                tfile_out<<std::setw(20)<<pops[ist]/SBZ;
            }

            indx=0;
            for(size_t ist=0; ist<Nstates; ist++){
                for(size_t jst=ist+1; jst<Nstates; jst++){
                    tfile_out<<std::setw(20)<<std::real(cohs[indx])/SBZ;
                    tfile_out<<std::setw(20)<<std::imag(cohs[indx])/SBZ;
                    indx++;
                }
            }

            ////////////////

            /// write conductivity

            ////////////////

            tfile_out << std::setw(20) << std::real(sigma[0]);
            tfile_out << std::setw(20) << std::imag(sigma[0]);
            tfile_out << std::setw(20) << std::abs(sigma[0]);

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
