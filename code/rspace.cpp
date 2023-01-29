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
typedef matrix<double> matrix_t;

#include <boost/format.hpp>

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];
        int tstep=std::stoi(argv[2])-1;

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

        ExternalField* E0=nullptr;

        ExternalField* Afield_x=new ExternalFieldFromData(Adata_x_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Afield_y=new ExternalFieldFromData(Adata_y_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Efield_x=new ExternalFieldFromData(Edata_x_fit,t0_fit,dt_fit,params.E0);
        ExternalField* Efield_y=new ExternalFieldFromData(Edata_y_fit,t0_fit,dt_fit,params.E0);

        auto tgrid=create_grid(params.tmin,params.tmax,params.Nt);
        double time=tgrid[tstep];

        //get vector potential at the given time step
        double Ax=(*Afield_x)(time);
        double Ay=(*Afield_y)(time);

        std::cout<<"Ax: "<<Ax<<" Ay: "<<Ay<<std::endl;

        HexagonalTBModel* tb=new HexagonalTBModel(params.a);

        //create kgrid
        Grid2D* kxygrid;
        Grid2D* Akxygrid;
        if(params.kgrid_type==kgrid_types::ucell){
            double Ox=0.;
            double Oy=0.;
            double b1x=2.*M_PI/(sqrt(3.)*params.a);
            double b1y=2.*M_PI/(params.a);
            double b2x=b1x;
            double b2y=-b1y;

            kxygrid=new UCellGrid2D(Ox,Oy,b1x,b1y,params.Nkx,b2x,b2y,params.Nky);
            Akxygrid=new UCellGrid2D(Ox+Ax,Oy+Ay,b1x,b1y,params.Nkx,b2x,b2y,params.Nky);
        }
        else{
            throw std::string("Integration is possible only for kgrid_type=ucell!");
        }

        Integrator2D* integrator_kxky=new Integrator2D(kxygrid);
        double SBZ=pow(2.*M_PI,2.)/(0.5*sqrt(3.)*params.a*params.a);

        //create graphene model
        Graphene* gm;
        if(params.model==models::hommelhoff){
            double e2p=params.e2p[0];
            double gamma=params.gamma[0];
            double s=params.s[0];
            gm=new GrapheneModel(params.a,e2p,gamma,s,params.Td,E0,E0);
        }
        else if(params.model==models::nlayer){
            gm=new NGraphene(tb,params.nlayers,
                        kxygrid,
                        params.e2p,params.gamma,params.s,params.Td,
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

        std::vector<matrix_t> dens_data(Nst);
        std::vector<matrix_t> coh_re_data(Nst*(Nst-1)/2);
        std::vector<matrix_t> coh_im_data(Nst*(Nst-1)/2);

        //read data from file
        std::string dens_t_fname=params.densfile_fname;
        replace(dens_t_fname,"%it",boost::str(boost::format("%06d") % (tstep+1)));

        std::cout<<dens_t_fname<<std::endl;

        size_t col=0;
        for(size_t ist=0; ist<Nst; ist++){
            dens_data[ist]=read_2D_from_file<matrix_t>(dens_t_fname,col,params.Nkx,params.Nky);
            col++;
        }

        size_t indx=0;
        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=ist+1; jst<Nst; jst++){
                coh_re_data[indx]=read_2D_from_file<matrix_t>(dens_t_fname,col,params.Nkx,params.Nky);
                col++;
                coh_im_data[indx]=read_2D_from_file<matrix_t>(dens_t_fname,col,params.Nkx,params.Nky);
                col++;
                indx++;
            }
        }

        //write rgrids to file
        /*std::ofstream rgrid_out(params.rgfile_fname);
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
        rgrid_out.close();*/

        MultiIndex indx_kxkyst({params.Nkx,params.Nky,Nst});
        size_t N_kxkyst=indx_kxkyst.size();

        MultiIndex indx_xyzst({params.Nx,params.Ny,params.Nz,Nst});
        size_t N_xyzst=indx_xyzst.size();

        //compute eigenstates in reciprocal space
        std::vector<vector_t> vecs(N_kxkyst);

        for(size_t ist=0; ist<Nst; ist++){
            for(size_t ikx=0; ikx<params.Nkx; ikx++){
                for(size_t iky=0; iky<params.Nky; iky++){
                    //we compute vectors using A-shifted grid!
                    double kx=(*Akxygrid)(ikx,iky)[0];
                    double ky=(*Akxygrid)(ikx,iky)[1];

                    size_t indx_ikxikyist=indx_kxkyst({ikx,iky,ist});
                    vecs[indx_ikxikyist]=gm->get_state(kx,ky,ist);
                }
            }
        }

        //create graphene material
        Orbital* pz=new Pzorb_normal(params.Z);
        Material graphene;

        //generate graphene
        for(size_t il=0; il<params.nlayers; il++){
            double x=0.;//shift of the layer depending on stacking
            
            stacking ABC=params.layers[il];
            switch(ABC){
                case stacking::A:
                    x=0.;
                    break;
                case stacking::B:
                    x=params.a/sqrt(3.);
                    break;
                case stacking::C:
                    x=2.*params.a/sqrt(3.);
                    break;
            }

            double z=il*params.d;//position of the layer in z coordinate

            //A..B atoms in graphene layer
            AtomsSet setA=GenerateGraphenePattern(pz,params.a,params.Nclx,params.Ncly,x,0.,z);
            AtomsSet setB=GenerateGraphenePattern(pz,params.a,params.Nclx,params.Ncly,x+params.a/sqrt(3.),0.,z);

            setA.compute_on_grid(kxygrid);
            setB.compute_on_grid(kxygrid);
            
            graphene.add_atomsset(setA);
            graphene.add_atomsset(setB);
        }

        //print atom positions to file
        /*std::ofstream atoms_out;
        atoms_out.open(params.afile_fname);
        graphene.print_atoms(atoms_out);
        atoms_out.close();*/

        MultiIndex indx_xyz({params.Nx,params.Ny,params.Nz});
        size_t N_xyz=indx_xyz.size();
        std::vector<double> Psi(N_xyz);

        //output real space data
        std::string rho_t_fname=params.rhofile_fname;
        replace(rho_t_fname,"%it",boost::str(boost::format("%06d") % (tstep+1)));
        std::ofstream rho_t_out(rho_t_fname);

        std::cout<<"Real space density will be written to: "<<rho_t_fname<<std::endl;

        rho_t_out<<std::scientific;
        rho_t_out<<std::setprecision(15);
        rho_t_out<<"#"<<std::setw(24)<<"rho_nocoh";
        rho_t_out<<     std::setw(25)<<"rho_coh";
        rho_t_out<<     std::setw(25)<<"rho_total";
        rho_t_out<<std::endl;   

        //3D real-space densities
        for(size_t ix=0; ix<params.Nx; ix++){
            std::cout<<"ix: "<<ix<<std::endl;
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<params.Nz; iz++){
                    double x=(*xygrid)(ix,iy)[0];
                    double y=(*xygrid)(ix,iy)[1];
                    double z=(*zgrid)[iz];

                    //integrate in reciprocal space                        
                    auto func=[x,y,z,
                        &dens_data,
                        &coh_re_data,&coh_im_data,
                        Nst,
                        &indx_kxkyst,
                        &vecs,&graphene](const size_t& ikx, const size_t& iky){

                        //we output:
                        //1..Nst A-modified eigenstates,
                        //rho_nocoh
                        vector<double> res(6);

                        //compute Bloch functions
                        std::vector<complex_t> BPhis(Nst);
                        std::vector<complex_t> dxBPhis(Nst);
                        std::vector<complex_t> dyBPhis(Nst);
                        for(size_t ist=0; ist<Nst; ist++){
                            BPhis[ist]=graphene.PhiI(ist,x,y,z,ikx,iky);
                            //dxBPhis[ist]=graphene.dPhiI(ist,x,y,z,ikx,iky,0);
                            //dyBPhis[ist]=graphene.dPhiI(ist,x,y,z,ikx,iky,1);
                        }

                        //compute eigenstates
                        std::vector<complex_t> Psis(Nst);
                        std::vector<complex_t> dxPsis(Nst);
                        for(size_t ist=0; ist<Nst; ist++){
                            complex_t Psi=0.;
                            complex_t dxPsi=0.;
                            for(size_t jst=0; jst<Nst; jst++){
                                size_t indx_ikxikyist=indx_kxkyst({ikx,iky,ist});
                                size_t indx_ikxikyjst=indx_kxkyst({ikx,iky,jst});
                                Psi+=vecs[indx_ikxikyist][jst]*BPhis[jst];
                                dxPsi+=vecs[indx_ikxikyist][jst]*dxBPhis[jst];
                            }
                            Psis[ist]=Psi;
                            dxPsis[ist]=dxPsi;
                        }

                        //compute non-coherent rho
                        double rho_nocoh=0.;
                        for(size_t ist=0; ist<Nst; ist++){
                            rho_nocoh+=dens_data[ist](ikx,iky)*std::norm(Psis[ist]);
                        }
                        res(0)=rho_nocoh;

                        //compute coherent rho
                        double rho_coh=0.;
                        size_t indx=0;
                        for(size_t ist=0; ist<Nst; ist++){
                            for(size_t jst=ist+1; jst<Nst; jst++){
                                double re=coh_re_data[indx](ikx,iky);
                                double im=coh_im_data[indx](ikx,iky);
                                complex_t coh=re+I*im;

                                complex_t psi_istjst=std::conj(Psis[ist])*Psis[jst];

                                rho_coh+=2.*std::real(coh*psi_istjst);

                                indx++;
                            }
                        }
                        res(1)=rho_coh;

                        //total rho
                        res(2)=rho_nocoh+rho_coh;

                        //compute non-coherent current
                        /*double Jx_nocoh=0.;
                        for(size_t ist=0; ist<Nst; ist++){
                            Jx_nocoh+=dens_data[ist](ikx,iky)*(
                                std::real(Psis[ist])*std::imag(dxPsis[ist])-std::imag(Psis[ist])*std::real(dxPsis[ist]));
                        }
                        res(3)=Jx_nocoh;*/

                        //compute coherent current
                        //res(4)=;
                        
                        //total current
                        //res(5)=;

                        return res;
                    };

                    vector<double> res=zero_vector<double>(6);
                    integrator_kxky->trapz(func,res);
                    res*=2./SBZ;//normalize for area of the BZ

                    size_t indx_ixiyiz=indx_xyz({ix,iy,iz});
                    Psi[indx_ixiyiz]=res(0);

                    //write data to file
                    for(size_t ist=0; ist<res.size(); ist++){
                        rho_t_out<<std::setw(25)<<res(ist);
                    }

                    rho_t_out<<std::endl;
                }
            }
        }
        rho_t_out.close();

        //integrate in real space
        //integrate in xy coordinates
        std::vector<double> res_xy(params.Nz);
        for(size_t iz=0; iz<zgrid->size(); iz++){
            auto int_xy=[iz,&Psi,&indx_xyz](const size_t& ix, const size_t& iy){
                size_t indx_ixiyiz=indx_xyz({ix,iy,iz});
                return Psi[indx_ixiyiz];
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

        std::cout<<"XYZkxky integral: "<<sum<<std::endl;

        /*for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<zgrid->size(); iz++){
                    double x=(*xygrid)(ix,iy)[0];
                    double y=(*xygrid)(ix,iy)[1];
                    double z=(*zgrid)[iz];

                    auto func=[x,y,z,
                    dens_vv,dens_cc,dens_cv_re,dens_cv_im,
                    &wfs,&wfs_g,
                    SBZ,
                    kxygrid](const size_t& ikx, const size_t& iky){
                        double kx=(*kxygrid)(ikx,iky)[0];
                        double ky=(*kxygrid)(ikx,iky)[1];

                        //complex_t psip=wfs.psip(x,y,z,kx,ky);
                        //complex_t psim=wfs.psim(x,y,z,kx,ky);

                        complex_t psip=wfs_g.psip(x,y,z,ikx,iky);
                        complex_t psim=wfs_g.psim(x,y,z,ikx,iky);

                        double rho_vv=std::norm(psip);
                        double rho_cc=std::norm(psim);
                        complex_t rho_vc=std::conj(psip)*psim;

                        complex_t dens_cv=dens_cv_re(ikx,iky)+I*dens_cv_im(ikx,iky);

                        double rho_nocoh=dens_vv(ikx,iky)*rho_vv+dens_cc(ikx,iky)*rho_cc;
                        double rho_coh=2.*std::real(dens_cv*rho_vc);
                        double rho_t=rho_nocoh+rho_coh;

                        //double rho_t=dens_vv(ikx,iky)*rho_vv
                        //            +dens_cc(ikx,iky)*rho_cc
                        //            +2.*std::real(dens_cv*rho_vc);

                        vector<double> res(5);
                        res(0)=rho_vv;
                        res(1)=rho_cc;
                        res(2)=rho_nocoh;
                        res(3)=rho_coh;
                        res(4)=rho_t;

                        //res(3)=(2./SBZ)*(rho_t-rho_vv);

                        //double res=(2./SBZ)*(rho_t-rho_vv);

                        //double res=rho_t;

                        //double res=(2./SBZ)*std::real(dens_cv*rho_vc);

                        //double res=(2./SBZ)*rho_t;

                        //double res=std::real(psim)-std::real(psim_g);

                        //double res=std::norm(psip)-std::norm(psim);
                        //double res=std::norm(psip_g)-std::norm(psim_g);

                        //double res=kx*kx+ky*ky;

                        //double res=pow(kx-ky,2.)*pow(cos(kx+ky),2.);

                        return res;
                    };

                    //double integr=kxygrid->integrate(func);
                    //std::cout<<integr<<std::endl;

                    //res(ix,iy,iz)(0)+=kxygrid->integrate(func);

                    integrator->trapz(func,res(ix,iy,iz));
                }
            }
        }
     
        for(size_t ix=0; ix<params.Nx; ix++){
            for(size_t iy=0; iy<params.Ny; iy++){
                for(size_t iz=0; iz<zgrid->size(); iz++){
                    for(size_t i=0; i<res(ix,iy,iz).size(); i++){
                        rho_t_out<<std::setw(25)<<res(ix,iy,iz)(i);
                    }
                    rho_t_out<<std::endl;
                }
            }
        }
        rho_t_out.close();

        /*double sum=0.;
        for(size_t iz=0; iz<zgrid->size(); iz++){
            std::function<double(const size_t& ix, const size_t& iy)> rho_xy=
            [res,iz](const size_t& ix, const size_t& iy){
                return res(ix,iy,iz);
            };
            sum+=xygrid->integrate(rho_xy);
        }
        std::cout<<"Sum: "<<sum<<std::endl;*/

    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}
