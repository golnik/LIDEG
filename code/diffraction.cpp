#include <iostream>
#include <memory>

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

#include "WFs.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;
typedef vector<double> vector_t;

#include <boost/format.hpp>

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

        ExternalField* Afield_x=new ExternalFieldFromData(Adata_x_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Afield_y=new ExternalFieldFromData(Adata_y_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Efield_x=new ExternalFieldFromData(Edata_x_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Efield_y=new ExternalFieldFromData(Edata_y_fit,t0_fit,dt_fit,params.E0);

        //diffraction spots
        typedef std::vector<std::vector<int>> BZ_t;
        BZ_t BZ0{{0,0}};
        BZ_t BZ1{{1,0},{1,1},{ 0,1},{-1, 0},{-1,-1},{0,-1}};
        BZ_t BZ2{{2,1},{1,2},{-1,1},{-2,-1},{-1,-2},{1,-1}};
        BZ_t BZ3{{2,0},{2,2},{ 0,2},{-2, 0},{-2,-2},{0,-2}};
        std::vector<BZ_t> zones{BZ0,BZ1,BZ2,BZ3};

        size_t nzones=zones.size();
        size_t nspots=0;
        for(size_t izone=0; izone<nzones; izone++){
            for(auto spot: zones[izone]){
                nspots++;
            }
        }

        //prepare grids
        auto tgrid=create_grid(params.tmin,params.tmax,params.Nt);

        //create kgrid
        double Ox=0.;
        double Oy=0.;
        double b1x=2.*M_PI/(sqrt(3.)*params.a);
        double b1y=2.*M_PI/(params.a);
        double b2x=b1x;
        double b2y=-b1y;

        vector_t O(2);
        O(0)=Ox;
        O(1)=Oy;

        vector_t a(2);
        a(0)=b1x;
        a(1)=b1y;

        vector_t b(2);
        b(0)=b2x;
        b(1)=b2y;

        Grid2D* kxygrid;
        if(params.kgrid_type==kgrid_types::ucell){
            kxygrid=new UCellGrid2D(Ox,Oy,b1x,b1y,params.Nkx,b2x,b2y,params.Nky);
        }
        else{
            throw std::string("Integration is possible only for kgrid_type=ucell!");
        }

        Integrator2D* integrator_kxky=new Integrator2D(kxygrid);
        double SBZ=pow(2.*M_PI,2.)/(0.5*sqrt(3.)*params.a*params.a);

        //create rgrid
        Grid2D* xygrid;
        if(params.rgrid_type==rgrid_types::ucell){
            double Ox=-(1./sqrt(3.))*params.a;
            double Oy=0.;
            double a1x=params.a/2.*sqrt(3.);
            double a1y=params.a/2.;
            double a2x=a1x;
            double a2y=-a1y;

            xygrid=new UCellGrid2D(Ox,Oy,a1x,a1y,params.Nx,a2x,a2y,params.Ny);
        }
        else{
            throw std::string("Integration is possible only for rgrid_type=ucell!");
        }

        Grid1D* zgrid=new RegularGrid1D(params.zmin,params.zmax,params.Nz);

        Integrator2D* integrator_xy=new Integrator2D(xygrid);
        Integrator1D* integrator_z=new Integrator1D(zgrid);

        //create required indices
        MultiIndex indx_xy({params.Nx,params.Ny});
        size_t N_xy=indx_xy.size();

        MultiIndex indx_xyz({params.Nx,params.Ny,params.Nz});
        size_t N_xyz=indx_xyz.size();

        //initialize output streams
        std::string diffr_intra_fname=params.diffr_fname;
        replace(diffr_intra_fname,"%type","intra");
        
        std::string diffr_inter_fname=params.diffr_fname;
        replace(diffr_inter_fname,"%type","inter");
        
        std::string diffr_total_fname=params.diffr_fname;
        replace(diffr_total_fname,"%type","total");

        std::vector<std::shared_ptr<std::ofstream>> ostreams{
            std::make_shared<std::ofstream>(diffr_intra_fname),
            std::make_shared<std::ofstream>(diffr_inter_fname),
            std::make_shared<std::ofstream>(diffr_total_fname)};

        //output header
        for(auto stream: ostreams){
            *stream<<std::scientific;
            *stream<<std::setprecision(8);

            *stream<<"#"<<std::setw(17)<<"Time[1]";

            *stream<<std::setw(18)<<"BZ0[2]"
                   <<std::setw(18)<<"BZ1[3]"
                   <<std::setw(18)<<"BZ2[4]"
                   <<std::setw(18)<<"BZ3[5]";

            size_t col_indx=6;
            for(auto BZ: zones){
                for(auto spot: BZ){
                    auto m=std::to_string(spot[0]);
                    auto n=std::to_string(spot[1]);
                    std::string label_str="BZ("+m+","+n+")["+std::to_string(col_indx)+"]";
                    *stream<<std::setw(18)<<label_str;
                    col_indx++;
                }
            }
            *stream<<std::endl;
        }

        //loop over time steps
        for(size_t it=0; it<params.Nt; it++){
            //read real space data
            std::string rho_t_fname=params.rhofile_fname;
            replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (it+1)));

            double time=tgrid[it];
            std::cout<<"Time step: "<<it+1<<std::endl;
            std::cout<<"Densities will be taken from "<<rho_t_fname<<" file"<<std::endl;

            //loop over intra, inter and total
            for(size_t icontrib=0; icontrib<3; icontrib++){
                //load data from file
                std::vector<double> dens_xyz;
                read_column_from_file(rho_t_fname,icontrib,dens_xyz);

                //create array for z-integrated data
                std::vector<double> dens_xy(N_xy);

                //we first integrate i z coordinate
                for(size_t ix=0; ix<params.Nx; ix++){
                    for(size_t iy=0; iy<params.Ny; iy++){
                        double res=0;
                        integrator_z->trapz(
                            [ix,iy,&dens_xyz,&indx_xyz](const size_t& iz){
                                size_t indx_ixiyiz=indx_xyz({ix,iy,iz});
                                return dens_xyz[indx_ixiyiz];
                            },res);
                        size_t indx_ixiy=indx_xy({ix,iy});
                        dens_xy[indx_ixiy]=res;
                    }
                }

                //data arrays for diffraction intensities
                std::vector<double> data(nspots);
                std::vector<double> avdata(nzones);

                //loop over Brillouin zones and spots in each zone
                size_t ispot=0;
                for(size_t izone=0; izone<nzones; izone++){
                    avdata[izone]=0.;
                    for(auto spot: zones[izone]){
                        auto m=spot[0];
                        auto n=spot[1];

                        //integrate in xy coordinates
                        auto int_xy=[m,n,O,a,b,&xygrid,dens_xy,&indx_xy,&params](const size_t& ix, const size_t& iy){
                            double x=(*xygrid)(ix,iy)[0];
                            double y=(*xygrid)(ix,iy)[1];
                            
                            size_t indx_ixiy=indx_xy({ix,iy});

                            double phi=0.;

                            //auto H=O+m*a+n*b;
                            //double H_r=H[0]*x+H[1]*y;
                            //std::complex<double> PW=exp(-I*(H_r-phi));

                            std::complex<double> PW=exp(I*2.*M_PI/params.a*(1./sqrt(3.)*(m+n)*x+(m-n)*y));

                            return dens_xy[indx_ixiy]*PW;
                        };

                        std::complex<double> res=0.;
                        integrator_xy->trapz(int_xy,res);

                        data[ispot]=std::abs(res);
                        avdata[izone]+=std::abs(res);

                        ispot++;
                    }
                }

                //output data to files
                *ostreams[icontrib]<<std::setw(18)<<time*au2fs;
                for(auto val: avdata){
                    *ostreams[icontrib]<<std::setw(18)<<val;
                }
                for(auto val: data){
                    *ostreams[icontrib]<<std::setw(18)<<val;
                }
                *ostreams[icontrib]<<std::endl;
            }
        }

        //close the streams
        for(auto stream: ostreams){
            stream->close();
        }
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}