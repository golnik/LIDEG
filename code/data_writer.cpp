#include <iostream>

#include "mini/ini.h"

#include "utils/utils.hpp"
#include "utils/grid.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "graphenemodel.hpp"
#include "model1/graphene.hpp"
#include "model1/WFs.hpp"
#include "model2/graphene2.hpp"
#include "Nlayer/nlayer.hpp"

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        ExternalField* E0=nullptr;

        HexagonalTBModel* tb=new HexagonalTBModel(params.a);

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

        //create graphene model
        Graphene* gm;
        if(params.model==models::hommelhoff){
            gm=new GrapheneModel(params.a,params.e2p,params.gamma,params.s,params.Td,E0,E0);
        }
        else if(params.model==models::nlayer){            
            double eps[]={params.e2p,params.e2p};
            double g[]={params.gamma,0.39/au2eV};
            double s[]={params.s,0.0};

            gm=new NGraphene(tb,params.nlayers,
                        kxygrid,
                        eps,g,s,
                        E0,E0);
        }

        //create rgrid
        double Ox=-(1./sqrt(3.))*params.a;
        double Oy=0.;
        double a1x=params.a/2.*sqrt(3.);
        double a1y=params.a/2.;
        double a2x=a1x;
        double a2y=-a1y;

        Grid2D* xygrid=new UCellGrid2D(Ox,Oy,a1x,a1y,params.Nx,a2x,a2y,params.Ny);
        Grid1D* zgrid=new RegularGrid1D(params.zmin,params.zmax,params.Nz);

        size_t Nst=gm->nstates();

        //write kgrid to file
        std::ofstream kgrid_out(params.kgfile_fname);
        kgrid_out<<params.Nkx*params.Nky<<std::endl;

        kgrid_out<<std::scientific;
        kgrid_out<<std::setprecision(8);

        for(size_t ikx=0; ikx<params.Nkx; ikx++){
            for(size_t iky=0; iky<params.Nky; iky++){
                kgrid_out<<std::setw(20)<<(*kxygrid)(ikx,iky)[0];
                kgrid_out<<std::setw(20)<<(*kxygrid)(ikx,iky)[1]<<std::endl;
            }
        }
        kgrid_out.close();

        //write rgrids to file
        std::ofstream rgrid_out(params.rgfile_fname);
        rgrid_out<<params.Nx*params.Ny*params.Nz<<std::endl;

        rgrid_out<<std::scientific;
        rgrid_out<<std::setprecision(8);

        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                rgrid_out<<std::setw(20)<<(*xygrid)(ix,iy)[0];
                rgrid_out<<std::setw(20)<<(*xygrid)(ix,iy)[1]<<std::endl;
            }
        }
        rgrid_out.close();

        //write kdata to file
        std::ofstream kout(params.pkfile_fname);

        size_t Nkx=kxygrid->size1();
        size_t Nky=kxygrid->size2();

        kout<<std::scientific;
        kout<<std::setprecision(8);

        kout<<"#";
        kout<<std::setw(19)<<"kx";
        kout<<std::setw(20)<<"ky";

        size_t col=3;
        for(size_t ist=0; ist<Nst; ist++){
            std::string Estr="E["+std::to_string(ist+1)+"]("+std::to_string(col)+")";
            kout<<std::setw(20)<<Estr;
            col++;
        }

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=ist+1; jst<Nst; jst++){
                std::string dx_str="dx["+std::to_string(ist+1)+","+std::to_string(jst+1)+"]("+std::to_string(col)+")";
                col++;
                std::string dy_str="dy["+std::to_string(ist+1)+","+std::to_string(jst+1)+"]("+std::to_string(col)+")";
                col++;

                kout<<std::setw(20)<<dx_str;
                kout<<std::setw(20)<<dy_str;
            }
        }
        kout<<std::endl;

        for(size_t ikx=0; ikx<Nkx; ikx++){
            for(size_t iky=0; iky<Nky; iky++){
                double kx=(*kxygrid)(ikx,iky)[0];
                double ky=(*kxygrid)(ikx,iky)[1];

                kout<<std::setw(20)<<kx<<std::setw(20)<<ky;

                for(size_t ist=0; ist<Nst; ist++){
                    kout<<std::setw(20)<<gm->get_energy(kx,ky,ist);
                }
                
                for(size_t ist=0; ist<Nst; ist++){
                    for(size_t jst=ist+1; jst<Nst; jst++){
                        double dip_x=gm->get_dipole(kx,ky,ist,jst,0);
                        double dip_y=gm->get_dipole(kx,ky,ist,jst,1);

                        kout<<std::setw(20)<<dip_x;
                        kout<<std::setw(20)<<dip_y;
                    }
                }

                kout<<std::endl;
            }
        }

        kout.close();

        //write rdata to file
        std::ofstream rout(params.prfile_fname);

        rout<<std::scientific;
        rout<<std::setprecision(8);

        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<params.Nz; iz++){
                    double x=(*xygrid)(ix,iy)[0];
                    double y=(*xygrid)(ix,iy)[1];
                    double z=(*zgrid)[iz];

                    rout<<std::setw(20)<<x;
                    rout<<std::setw(20)<<y;
                    rout<<std::setw(20)<<z;

                    rout<<std::endl;
                }
            }
        }

        rout.close();
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}