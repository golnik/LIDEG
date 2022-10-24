#include <iostream>
#include <complex>

#include "utils/utils.hpp"
#include "external_field.hpp"
#include "model1/graphene.hpp"
#include "Nlayer/nlayer.hpp"

int main(int argc, char** argv){
    double a=2.46/au2A;
    HexagonalTBModel* tb=new HexagonalTBModel(a);

    double e2p=0.0;
    double gamma=-3.033/au2eV;
    double s0=0.129;

    double eps[]={e2p,e2p};
    double g[]={gamma,0.39/au2eV};
    double s[]={s0,0.0};
    size_t Nlayers=2;

    ExternalField* E0=nullptr;

    //create graphene model
    GrapheneModel gm(a,e2p,gamma,s0,E0,E0);

    NGraphene ngr(tb,Nlayers,eps,g,s,E0,E0);

    size_t Nst=2*Nlayers;

    double kymin=-1.4;
    double kymax=1.4;
    size_t Nky=10001;

    double kx=0.0;

    std::cout<<"# ky ";
    for(size_t ist=0; ist<Nst; ist++){
        for(size_t jst=ist+1; jst<Nst; jst++){
            std::cout<<"<"<<ist<<"|"<<jst<<">"<<" ";
        }
    }
    std::cout<<std::endl;

    double dky=(kymax-kymin)/(Nky-1);
    double ky=kymin;
    for(size_t iky=0; iky<Nky; iky++){
        ngr.solve(kx,ky);

        auto evals=ngr.getEnergies();

        std::cout<<ky<<" ";
        
        for(size_t ist=0; ist<Nst; ist++){
            for(size_t jst=ist+1; jst<Nst; jst++){
                std::cout<<std::abs(ngr.getDipole(ist,jst))<<" ";
            }
        }

        std::cout<<gm.dx(kx,ky)<<" ";

        /*for(size_t ist=0; ist<Nst; ist++){
            std::cout<<au2eV*evals[ist]<<" ";
        }*/

        //std::cout<<gm.ep(kx,ky)*au2eV<<" "<<gm.em(kx,ky)*au2eV<<" ";
        //std::cout<<evals(0)<<" "<<evals(1)<<" ";

        //std::cout<<evals<<std::endl;
        //std::cout<<evecs<<std::endl;//.transpose()*evecs(1)<<std::endl;
        //std::cout<<evecs.col(0)<<std::endl;

        //auto e0=evals(0);
        //auto e1=evals(1);
        //auto vec0=evecs.col(0);
        //auto vec1=evecs.col(1);

        //std::cout<<"norms: "<<std::endl;

        //std::cout<<vec0.dot(vec0)<<std::endl;
        //std::cout<<vec1.dot(vec0)<<std::endl;

        //complex_t dip=I*vec0.dot(S.inverse()*dH0*vec1)/(e0-e1);

        //std::cout<<std::real(dip)<<" "<<std::imag(dip)<<" ";
        //std::cout<<gm.dy(kx,ky)<<" ";

        //std::cout<<evecs.adjoint().col(0).dot(evecs.col(1))<<std::endl;

        //std::cout<<std::real(fval)<<" "<<std::imag(fval)<<" ";
        //std::cout<<std::real(dfval)<<" "<<std::imag(dfval)<<" ";
        
        std::cout<<std::endl;

        ky+=dky;
    }

    return 0;
}