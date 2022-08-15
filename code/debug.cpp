#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;

#include "mini/ini.h"

#include "utils/utils.hpp"
#include "parser.hpp"
#include "external_field.hpp"
#include "model1/graphene.hpp"
#include "model1/WFs.hpp"
#include "model2/graphene2.hpp"

/*#include <boost/math/differentiation/autodiff.hpp>
using namespace boost::math::differentiation;
using namespace boost::math::differentiation::detail;

class NewField:
public ExternalFieldFromData{
public:
    NewField(const ExternalFieldFromData& ef):
    ExternalFieldFromData(ef){}

    template<typename RealType, size_t Order>
    fvar<RealType,Order> operator()(fvar<RealType, Order> const& cr) const{
        using root_type=typename fvar<RealType,Order>::root_type;
        constexpr size_t order=fvar<RealType,Order>::order_sum;

        //zero order derivative
        root_type const d0=ExternalFieldFromData::operator()(static_cast<root_type>(cr));

        if constexpr (order==0){
            return fvar<RealType,Order>(d0);
        }
        else{
            root_type derivatives[order]{0.0};
            derivatives[0]=d0;
            derivatives[1]=ExternalFieldFromData::derivative(static_cast<root_type>(cr));
            return cr.apply_derivatives(order,[&derivatives](size_t i){return derivatives[i];});
        }
    }

    double operator()(const double& t) const{
        return NewField::operator()(make_fvar<double,0>(t)).derivative(0);
        //return 0.;
    }
};

class GrapheneModel2{
    typedef std::vector<complex_t> state_type;    
public:
    GrapheneModel2(const double& a,
                   const double& e2p, const double& gamma, const double& s,
                   ExternalField* Ex, ExternalField* Ey,
                   ExternalField* Ax, ExternalField* Ay):
    _a(a),
    _e2p(e2p),_g0(gamma),_s0(s){
        _Ex=new NewField(*static_cast<ExternalFieldFromData*>(Ex));
        _Ey=new NewField(*static_cast<ExternalFieldFromData*>(Ey));
        _Ax=new NewField(*static_cast<ExternalFieldFromData*>(Ax));
        _Ay=new NewField(*static_cast<ExternalFieldFromData*>(Ay));
    }

    template<typename T>
    T h1(const T& kx, const T& ky) const{
        return 2.*cos(0.5*ky*_a/sqrt(3.))*cos(0.5*kx*_a);
    }

    template<typename T>
    T h2(const T& kx, const T& ky) const{
        return sin(ky*_a/sqrt(3.))-2.*sin(0.5*ky*_a/sqrt(3.))*cos(0.5*kx*_a);
    }

    template<typename T>
    T chi(const T& kx, const T& ky) const{
        return atan(h2(kx,ky)/h1(kx,ky));
    }

    template<typename T1, typename T2, typename T3>
    promote<T1,T2,T3> chi_t(const T1& kx0, const T2& ky0, const T3& t) const{
        auto kxt=kx0+(*_Ax)(t);
        auto kyt=ky0+(*_Ay)(t);
        return chi(kxt,kyt);
    }

    double dchidt(const double& kx0, const double& ky0, const double& t) const{
        auto const S=make_fvar<double,1>(t);
        auto const autodiff=chi_t(kx0,ky0,S);
        return autodiff.derivative(1);
    }

    double f(const double& kx, const double& ky) const{
        return sqrt(h1(kx,ky)*h1(kx,ky)+h2(kx,ky)*h2(kx,ky));
    }

    double E1(const double& kx, const double& ky) const{
        return (_e2p-_g0*f(kx,ky))/(1.+_s0*f(kx,ky));
    }

    double E2(const double& kx, const double& ky) const{
        return (_e2p+_g0*f(kx,ky))/(1.-_s0*f(kx,ky));
    }

    void propagate(const state_type &c, state_type &dcdt, const double t,
                   const double& kx0, const double& ky0) const{
        double kxt=kx0+(*_Ax)(t);
        double kyt=ky0+(*_Ay)(t);

        double dchi=dchidt(kx0,ky0,t);

        double V11=(_e2p-_g0*f(kxt,kyt))/(1.+_s0*f(kxt,kyt)-0.5*dchi);
        double V12=0.5*dchi;
        double V21=0.5*dchi;
        double V22=(_e2p+_g0*f(kxt,kyt))/(1.-_s0*f(kxt,kyt)-0.5*dchi);

        dcdt[0]=-I*(V11*c[0]+V12*c[1]);
        dcdt[1]=-I*(V21*c[0]+V22*c[1]);

        return;
    }
private:
    double _e2p = 0.0;
    double _g0  = 3.033/au2eV;
    double _s0  = 0.129;
    double _a   = 0.246/au2nm;

    NewField* _Ex;
    NewField* _Ey;
    NewField* _Ax;
    NewField* _Ay;
};*/

#include <boost/math/quadrature/gauss.hpp>
using namespace boost::math::quadrature;

#define __cplusplus
#include "utils/gauss_legendre.h"


/* n = 16 */
//static double x16[8] = {0.0950125098376374401853193,0.2816035507792589132304605,0.4580167776572273863424194,0.6178762444026437484466718,0.7554044083550030338951012,0.8656312023878317438804679,0.9445750230732325760779884,0.9894009349916499325961542};
//static double w16[8] = {0.1894506104550684962853967,0.1826034150449235888667637,0.1691565193950025381893121,0.1495959888165767320815017,0.1246289712555338720524763,0.0951585116824927848099251,0.0622535239386478928628438,0.0271524594117540948517806};

/* n = 128 */
//static double x128[64] = {0.0122236989606157641980521,0.0366637909687334933302153,0.0610819696041395681037870,0.0854636405045154986364980,0.1097942311276437466729747,0.1340591994611877851175753,0.1582440427142249339974755,0.1823343059853371824103826,0.2063155909020792171540580,0.2301735642266599864109866,0.2538939664226943208556180,0.2774626201779044028062316,0.3008654388776772026671541,0.3240884350244133751832523,0.3471177285976355084261628,0.3699395553498590266165917,0.3925402750332674427356482,0.4149063795522750154922739,0.4370245010371041629370429,0.4588814198335521954490891,0.4804640724041720258582757,0.5017595591361444642896063,0.5227551520511754784539479,0.5434383024128103634441936,0.5637966482266180839144308,0.5838180216287630895500389,0.6034904561585486242035732,0.6228021939105849107615396,0.6417416925623075571535249,0.6602976322726460521059468,0.6784589224477192593677557,0.6962147083695143323850866,0.7135543776835874133438599,0.7304675667419088064717369,0.7469441667970619811698824,0.7629743300440947227797691,0.7785484755064119668504941,0.7936572947621932902433329,0.8082917575079136601196422,0.8224431169556438424645942,0.8361029150609068471168753,0.8492629875779689691636001,0.8619154689395484605906323,0.8740527969580317986954180,0.8856677173453972174082924,0.8967532880491581843864474,0.9073028834017568139214859,0.9173101980809605370364836,0.9267692508789478433346245,0.9356743882779163757831268,0.9440202878302201821211114,0.9518019613412643862177963,0.9590147578536999280989185,0.9656543664319652686458290,0.9717168187471365809043384,0.9771984914639073871653744,0.9820961084357185360247656,0.9864067427245862088712355,0.9901278184917343833379303,0.9932571129002129353034372,0.9957927585349811868641612,0.9977332486255140198821574,0.9990774599773758950119878,0.9998248879471319144736081};
//static double w128[64] = {0.0244461801962625182113259,0.0244315690978500450548486,0.0244023556338495820932980,0.0243585572646906258532685,0.0243002001679718653234426,0.0242273192228152481200933,0.0241399579890192849977167,0.0240381686810240526375873,0.0239220121367034556724504,0.0237915577810034006387807,0.0236468835844476151436514,0.0234880760165359131530253,0.0233152299940627601224157,0.0231284488243870278792979,0.0229278441436868469204110,0.0227135358502364613097126,0.0224856520327449668718246,0.0222443288937997651046291,0.0219897106684604914341221,0.0217219495380520753752610,0.0214412055392084601371119,0.0211476464682213485370195,0.0208414477807511491135839,0.0205227924869600694322850,0.0201918710421300411806732,0.0198488812328308622199444,0.0194940280587066028230219,0.0191275236099509454865185,0.0187495869405447086509195,0.0183604439373313432212893,0.0179603271850086859401969,0.0175494758271177046487069,0.0171281354231113768306810,0.0166965578015892045890915,0.0162550009097851870516575,0.0158037286593993468589656,0.0153430107688651440859909,0.0148731226021473142523855,0.0143943450041668461768239,0.0139069641329519852442880,0.0134112712886163323144890,0.0129075627392673472204428,0.0123961395439509229688217,0.0118773073727402795758911,0.0113513763240804166932817,0.0108186607395030762476596,0.0102794790158321571332153,0.0097341534150068058635483,0.0091830098716608743344787,0.0086263777986167497049788,0.0080645898904860579729286,0.0074979819256347286876720,0.0069268925668988135634267,0.0063516631617071887872143,0.0057726375428656985893346,0.0051901618326763302050708,0.0046045842567029551182905,0.0040162549837386423131943,0.0034255260409102157743378,0.0028327514714579910952857,0.0022382884309626187436221,0.0016425030186690295387909,0.0010458126793403487793129,0.0004493809602920903763943};

double f(double x, double y, void* data){
    return (1.-x*x-y*y)*exp(-x*x);
}

double f(double x, void* data){
	return sin(x);
}

std::vector<double> create_quadrature_grid(
    const double& a, const double& b, const size_t& N){
    std::vector<double> grid(N);

    double A=0.5*(b-a);
    double B=0.5*(b+a);

    for(size_t i=0; i<N/2; i++){
        double xf=B+A*x128[i];
        double xb=B-A*x128[i];

        grid[N/2+i]  =xf;
        grid[N/2-i-1]=xb;
    }

    return grid;
}

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

typedef matrix<double> matrix_t;

int main(int argc, char** argv){
    try{
        std::string fname=argv[1];

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

        double a=-3.;
        double b=4.;
        double c=-4.;
        double d=5.;
        size_t N=128;

        std::cout<<std::fixed;
        std::cout<<std::setprecision(6);

        std::vector<double> xgrid=create_quadrature_grid(a,b,N);
        std::vector<double> ygrid=create_quadrature_grid(c,d,N);

        matrix_t fgrid(N,N);

        for(size_t ix=0; ix<N; ix++){
            for(size_t iy=0; iy<N; iy++){
                fgrid(ix,iy)=f(xgrid[ix],ygrid[iy],NULL);
            }
        }

        double res=0.;
        for(size_t ix=0; ix<N/2; ix++){
            size_t ixf=N/2+ix;
            size_t ixb=N/2-ix-1;

            for(size_t iy=0; iy<N/2; iy++){
                size_t iyf=N/2+iy;
                size_t iyb=N/2-iy-1;

                res+=w128[ix]*w128[iy]*(
                    fgrid(ixf,iyf)+fgrid(ixf,iyb)
                   +fgrid(ixb,iyf)+fgrid(ixb,iyb)
                );
            }
        }
        res*=0.5*(b-a)*0.5*(d-c);

        std::cout<<std::endl;
        
        std::cout<<res<<std::endl;

        double approx=gauss_legendre_2D_cube(N,f,NULL,a,b,c,d);
        //double approx=gauss_legendre(N,f,NULL,a,b);

        std::cout<<approx<<std::endl;

        //read fields from file
        /*std::vector<double> tgrid_fit;
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

        double kmin=-20.*au2nm;
        double kmax= 20.*au2nm;

        auto kx_grid=create_grid(kmin,kmax,params.Nkx);
        auto ky_grid=create_grid(kmin,kmax,params.Nky);

        fs::create_directories("./debug");
        std::ofstream debug_out("debug/debug.out");

        for(size_t ikx=0; ikx<params.Nkx; ikx++){
            for(size_t iky=0; iky<params.Nky; iky++){
                double kx0=kx_grid[ikx];
                double ky0=ky_grid[iky];

                //debug_out<<" ";
                
                debug_out<<std::abs(gm.f(kx0,ky0))<<" ";
                debug_out<<std::abs(gm2.f(kx0,ky0))<<" ";
                
                debug_out<<std::endl;
            }
        }
        
        debug_out.close();*/
    }catch(std::string er){
        std::cout<<' '<<er<<std::endl;
        std::cout<<" Task not accomplished.\n";
        return 1;
    }
    std::cout<<"\n Tasks accomplished.\n";
    return 0;
}