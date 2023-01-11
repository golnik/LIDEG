#include <iostream>

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
typedef vector<complex_t> vector_t;

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        //create output directory if does not exist
        fs::path outpath(params.outdir);
        fs::create_directory(outpath);

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

        Integrator2D* integrator_kxky=new Integrator2D(kxygrid);
        double SBZ=pow(2.*M_PI,2.)/(0.5*sqrt(3.)*params.a*params.a);

        //create graphene model
        Graphene* gm;
        if(params.model==models::hommelhoff){
            gm=new GrapheneModel(params.a,params.e2p,params.gamma,params.s,params.Td,E0,E0);
        }
        else if(params.model==models::nlayer){            
            std::vector<double> eps={params.e2p,params.e2p};
            std::vector<double> g={params.gamma,0.39/au2eV};
            std::vector<double> s={params.s,0.0};

            gm=new NGraphene(tb,params.nlayers,
                        kxygrid,
                        eps,g,s,
                        E0,E0);
        }

        //create rgrid
        Grid2D* xygrid;
        if(params.rgrid_type==rgrid_types::rectan){
            Grid1D* xgrid=new RegularGrid1D(params.xmin,params.xmax,params.Nx);
            Grid1D* ygrid=new RegularGrid1D(params.ymin,params.ymax,params.Ny);

            xygrid=new RegularGrid2D(xgrid,ygrid);
        }
        else if(params.rgrid_type==rgrid_types::ucell){
            double Ox=-(1./sqrt(3.))*params.a;
            double Oy=0.;
            double a1x=params.a/2.*sqrt(3.);
            double a1y=params.a/2.;
            double a2x=a1x;
            double a2y=-a1y;

            xygrid=new UCellGrid2D(Ox,Oy,a1x,a1y,params.Nx,a2x,a2y,params.Ny);
        }

        Grid1D* zgrid=new RegularGrid1D(params.zmin,params.zmax,params.Nz);

        Integrator2D* integrator_xy=new Integrator2D(xygrid);
        Integrator1D* integrator_z=new Integrator1D(zgrid);

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
                for(size_t iz=0; iz<params.Nz; iz++){
                    double x=(*xygrid)(ix,iy)[0];
                    double y=(*xygrid)(ix,iy)[1];
                    double z=(*zgrid)[iz];

                    rgrid_out<<std::setw(20)<<x;
                    rgrid_out<<std::setw(20)<<y;
                    rgrid_out<<std::setw(20)<<z;
                    rgrid_out<<std::endl;
                }
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
                    kout<<std::setw(20)<<gm->get_energy_grad(kx,ky,ist,0);
                    kout<<std::setw(20)<<gm->get_energy_grad(kx,ky,ist,1);
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

        MultiIndex indx_kxkyst({params.Nkx,params.Nky,Nst});
        size_t N_kxkyst=indx_kxkyst.size();

        MultiIndex indx_xyzst({params.Nx,params.Ny,params.Nz,Nst});
        size_t N_xyzst=indx_xyzst.size();

        //compute eigenstates in reciprocal space
        std::vector<vector_t> vecs(N_kxkyst);

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t ikx=0; ikx<Nkx; ikx++){
                for(size_t iky=0; iky<Nky; iky++){
                    double kx=(*kxygrid)(ikx,iky)[0];
                    double ky=(*kxygrid)(ikx,iky)[1];

                    size_t indx_ikxikyist=indx_kxkyst({ikx,iky,ist});
                    vecs[indx_ikxikyist]=gm->get_state(kx,ky,ist);
                }
            }
        }

        //create graphene material
        Orbital* pz=new Pzorb_normal(params.Z);

        //double l=3.46/au2A;

        AtomsSet setA=GenerateGraphenePattern(pz,params.a,params.Nclx,params.Ncly,0.,0.,0.);
        AtomsSet setB=GenerateGraphenePattern(pz,params.a,params.Nclx,params.Ncly,params.a/sqrt(3.),0.,0.);

        setA.compute_on_grid(kxygrid);
        setB.compute_on_grid(kxygrid);

        //AtomsSet setAA=GenerateGraphenePattern(pz,params.a,params.Nclx,params.Ncly,params.a/sqrt(3.),0.,l);
        //AtomsSet setBB=GenerateGraphenePattern(pz,params.a,params.Nclx,params.Ncly,0.,0.,l);


        Material graphene;
        graphene.add_atomsset(setA);
        graphene.add_atomsset(setB);

        //graphene.add_atomsset(setAA);
        //graphene.add_atomsset(setBB);

        //print atom positions to file
        std::ofstream atoms_out;
        atoms_out.open(params.afile_fname);
        graphene.print_atoms(atoms_out);
        atoms_out.close();

        std::vector<double> Psi(N_xyzst);//eigenstates data in real space

        //write rdata to file
        std::ofstream rout(params.prfile_fname);
        rout<<std::scientific;
        rout<<std::setprecision(15);

        for(size_t ix=0; ix<params.Nx; ix++){
            std::cout<<"ix: "<<ix<<std::endl;
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<params.Nz; iz++){
                    double x=(*xygrid)(ix,iy)[0];
                    double y=(*xygrid)(ix,iy)[1];
                    double z=(*zgrid)[iz];

                    //integrate in reciprocal space                        
                    auto func=[x,y,z,
                        Nst,
                        &indx_kxkyst,
                        &vecs,&graphene](const size_t& ikx, const size_t& iky){

                        //we will return eigenstates
                        vector<double> res(Nst);

                        //compute Bloch functions
                        std::vector<complex_t> BPhis(Nst);
                        for(size_t ist=0; ist<Nst; ist++){
                            BPhis[ist]=graphene.PhiI(ist,x,y,z,ikx,iky);
                        }

                        //compute eigenstates
                        for(size_t ist=0; ist<Nst; ist++){
                            complex_t Psi=0.;
                            for(size_t jst=0; jst<Nst; jst++){
                                size_t indx_ikxikyist=indx_kxkyst({ikx,iky,ist});
                                size_t indx_ikxikyjst=indx_kxkyst({ikx,iky,jst});
                                Psi+=vecs[indx_ikxikyist][jst]*BPhis[jst];
                            }
                            res(ist)=std::norm(Psi);
                        }
                        
                        return res;
                    };

                    vector<double> res=zero_vector<double>(Nst);
                    integrator_kxky->trapz(func,res);
                    res*=2./SBZ;//normalize for area of the BZ

                    for(size_t ist=0; ist<Nst; ist++){
                        size_t indx_ixiyizist=indx_xyzst({ix,iy,iz,ist});
                        Psi[indx_ixiyizist]=res(ist);
                    }

                    //write eigenstates to files
                    for(size_t ist=0; ist<res.size(); ist++){
                        rout<<std::setw(25)<<res(ist);
                    }

                    rout<<std::endl;
                }
            }
        }
        rout.close();

        //integrate in real space
        for(size_t ist=0; ist<Nst; ist++){
            //integrate in xy coordinates
            std::vector<double> res_xy(params.Nz);
            for(size_t iz=0; iz<zgrid->size(); iz++){
                auto int_xy=[iz,ist,&Psi,&indx_xyzst](const size_t& ix, const size_t& iy){
                    size_t indx_ixiyizist=indx_xyzst({ix,iy,iz,ist});
                    return Psi[indx_ixiyizist];
                };
                double res=0.;
                integrator_xy->trapz(int_xy,res);
                res_xy[iz]=res;
            }

            //integrate in z coordinate
            double sum=0;
            integrator_z->trapz(
                [&res_xy](const size_t& iz){
                    return res_xy[iz];
                },sum);

            std::cout<<"XYZkxky integral of eigenstate "<<ist+1<<": "<<sum<<std::endl;
        }

    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}