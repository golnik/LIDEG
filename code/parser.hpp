#ifndef PARSER_HPP
#define PARSER_HPP

#include "mini/ini.h"
//#include "utils/utils.hpp"

#include <string>
#include <fstream>
#include <regex>
#include <iomanip>

#include <filesystem>
namespace fs=std::filesystem;

enum class kgrid_types: unsigned int {ucell=0, kpoint};
enum class rgrid_types: unsigned int {ucell=0, rectan};
enum class models: unsigned int {hommelhoff=0, nlayer};
enum class stacking: unsigned int {A=0, B, C};

struct Parameters{
    double tmin;
    double tmax;
    size_t Nt;

    double ddt;
    double err_abs;
    double err_rel;

    double dkx;
    size_t Nkx;

    double dky;
    size_t Nky;

    kgrid_types kgrid_type;

    size_t Nclx;
    size_t Ncly;
    double Rmax;

    double xmin;
    double xmax;
    size_t Nx;

    double ymin;
    double ymax;
    size_t Ny;

    double zmin;
    double zmax;
    size_t Nz;        

    rgrid_types rgrid_type;

    double E0;
    std::string field_fname;

    models model;
    size_t nlayers;
    std::vector<stacking> layers;

    double d;
    double a;
    std::vector<double> e2p;
    std::vector<double> gamma;
    std::vector<double> s;
    double Z;
    double Td;

    std::string outdir;
    std::string kgfile_fname;
    std::string rgfile_fname;
    std::string tfile_fname;
    std::string densfile_fname;
    std::string rhofile_fname;
    std::string afile_fname;
    std::string pkfile_fname;
    std::string prfile_fname;

    void print(std::ostream& out) const{
        out<<"********************************************\n";
        out<<"********* PROPAGATION PARAMETERS: **********\n";
        out<<"********************************************\n";
        out<<std::fixed;
        out<<std::setprecision(5);
        out<<"\ttmin: "<<tmin*au2fs<<", fs\n";
        out<<"\ttmax: "<<tmax*au2fs<<", fs\n";
        out<<"\tNt: "<<Nt<<std::endl;
        out<<std::scientific;
        out<<"\tddt: "<<ddt*au2fs<<", fs\n";
        out<<"\tabsolute error: "<<err_abs<<std::endl;
        out<<"\trelative error: "<<err_rel<<std::endl;

        out<<"********************************************\n";
        out<<"************ FIELD PARAMETERS: *************\n";
        out<<"********************************************\n";
        out<<"\tfield will be loaded from: "<<field_fname<<std::endl;
        out<<std::scientific;
        out<<"\tfield intensity: "<<E0*au2Vnm<<", V/nm\n";

        out<<"********************************************\n";
        out<<"********* PARAMETERS OF THE MODEL: *********\n";
        out<<"********************************************\n";
        out<<std::fixed;
        out<<std::setprecision(5);

        out<<"\tmodel name: ";
        switch(model){
            case models::hommelhoff:
                out<<"Hommelhoff model\n";
                break;
            case models::nlayer:
                out<<"nlayer numerical model\n";
                break;
        }

        out<<"\tlattice constant: "<<a*au2A<<", A"<<std::endl;
        out<<"\tnumber of layers: "<<nlayers<<std::endl;

        if(nlayers!=1){
            out<<"\tinterlayer distance: "<<d*au2A<<", A"<<std::endl;
            out<<"\tlayers stacking: ";
            for(auto it: layers){
                switch(it){
                    case stacking::A:
                        out<<"A,";
                        break;
                    case stacking::B:
                        out<<"B,";
                        break;
                    case stacking::C:
                        out<<"C,";
                        break;
                }
            }
            out<<"\n";
        }

        out<<"\torbital energy: ";
        for(auto it: e2p){
            out<<it*au2eV<<", ";
        }
        out<<"eV\n";

        out<<"\thopping energy: ";
        for(auto it: gamma){
            out<<it*au2eV<<", ";
        }
        out<<"eV\n";

        out<<"\toverlap integral: ";
        for(auto it: s){
            out<<it<<", ";
        }
        out<<"\n";

        out<<"\torbital charge: "<<Z<<std::endl;
        out<<"\tcoherence time: "<<Td*au2fs<<", fs\n";

        return;
    }
};

class Parser{
public:
    Parser(Parameters& params):
    _params(params){}
    ~Parser(){}

    void analyze(const std::string& fname){
        mINI::INIFile file(fname);
        mINI::INIStructure ini;

        file.read(ini);

        //parse propagator section
        std::string tmin_str=ini.get("propagator").get("tmin");
        _params.tmin=std::stod(tmin_str)/au2fs;

        std::string tmax_str=ini.get("propagator").get("tmax");
        _params.tmax=std::stod(tmax_str)/au2fs;

        std::string Nt_str=ini.get("propagator").get("Nt");
        _params.Nt=std::stoi(Nt_str);

        std::string ddt_str=ini.get("propagator").get("ddt");
        _params.ddt=std::stod(ddt_str)/au2fs;

        std::string err_abs_str=ini.get("propagator").get("err_abs");
        _params.err_abs=std::stod(err_abs_str);

        std::string err_rel_str=ini.get("propagator").get("err_rel");
        _params.err_rel=std::stod(err_rel_str);

        //parse kgrid section
        std::string dkx_str=ini.get("kgrid").get("dkx");
        _params.dkx=std::stod(dkx_str)*au2nm;

        std::string Nkx_str=ini.get("kgrid").get("Nkx");
        _params.Nkx=std::stoi(Nkx_str);  

        std::string dky_str=ini.get("kgrid").get("dky");
        _params.dky=std::stod(dky_str)*au2nm;

        std::string Nky_str=ini.get("kgrid").get("Nky");
        _params.Nky=std::stoi(Nky_str);

        std::string kgrid_type_str=ini.get("kgrid").get("type");
        trim(kgrid_type_str);
        if(kgrid_type_str=="kpoint"){
            _params.kgrid_type=kgrid_types::kpoint;
        }
        else if(kgrid_type_str=="ucell"){
            _params.kgrid_type=kgrid_types::ucell;
        }
        else{
            throw std::string("The requested kgrid_type does not exist!");
        }

        //parse rgrid section
        std::string Nclx_str=ini.get("rgrid").get("Nclx");
        _params.Nclx=std::stoi(Nclx_str);
        std::string Ncly_str=ini.get("rgrid").get("Ncly");
        _params.Ncly=std::stoi(Ncly_str);
        std::string Rmax_str=ini.get("rgrid").get("Rmax");
        _params.Rmax=std::stod(Rmax_str)/au2A;

        _params.xmin=std::stod(ini.get("rgrid").get("xmin"))/au2A;
        _params.xmax=std::stod(ini.get("rgrid").get("xmax"))/au2A;
        _params.Nx=std::stoi(ini.get("rgrid").get("Nx"));

        _params.ymin=std::stod(ini.get("rgrid").get("ymin"))/au2A;
        _params.ymax=std::stod(ini.get("rgrid").get("ymax"))/au2A;
        _params.Ny=std::stoi(ini.get("rgrid").get("Ny"));

        _params.zmin=std::stod(ini.get("rgrid").get("zmin"))/au2A;
        _params.zmax=std::stod(ini.get("rgrid").get("zmax"))/au2A;
        _params.Nz=std::stoi(ini.get("rgrid").get("Nz"));

        std::string rgrid_type_str=ini.get("rgrid").get("type");
        trim(rgrid_type_str);
        if(rgrid_type_str=="rectan"){
            _params.rgrid_type=rgrid_types::rectan;
        }
        else if(rgrid_type_str=="ucell"){
            _params.rgrid_type=rgrid_types::ucell;
        }
        else{
            throw std::string("The requested rgrid_type does not exist!");
        }

        //parse field file
        _params.E0=std::stod(ini.get("field").get("E0"))/au2Vnm;
        _params.field_fname=ini.get("field").get("fname");

        //parse model parameters
        std::string model_str=ini.get("system").get("model");
        trim(model_str);
        if(model_str=="hommelhoff"){
            _params.model=models::hommelhoff;
        }
        else if(model_str=="nlayer"){
            _params.model=models::nlayer;
        }
        else{
            throw std::string("The requested model does not exist!");
        }

        std::string layers_str=ini.get("system").get("layers");
        if(!layers_str.empty()){
            std::regex pattern("\\w+");
            for(auto it=std::sregex_iterator(layers_str.begin(),layers_str.end(),pattern);
                     it!=std::sregex_iterator(); it++){
                std::smatch m=*it;
                std::string lstr=m.str();
                if(lstr=="A"){
                    _params.layers.push_back(stacking::A);
                }
                else if(lstr=="B"){
                    _params.layers.push_back(stacking::B);
                }
                else if(lstr=="C"){
                    _params.layers.push_back(stacking::C);
                }
                else{
                    throw std::string("The requested layers stacking does not exist!");
                }
            }    
        }
        else{
            _params.layers.push_back(stacking::A);
        }
        _params.nlayers=_params.layers.size();

        std::string d_str=ini.get("system").get("d");
        _params.d=std::stod(d_str)/au2A;

        std::string a_str=ini.get("system").get("a");
        _params.a=std::stod(a_str)/au2A;

        std::string e2p_str=ini.get("system").get("e2p");
        if(!e2p_str.empty()){
            _params.e2p=parse_array<double>(e2p_str);
            for(size_t i=0; i<_params.e2p.size(); i++){
                _params.e2p[i]/=au2eV;
            }
        }
        else{
            throw std::string("e2p parameter is not specified!");
        }

        std::string gamma_str=ini.get("system").get("gamma");
        if(!gamma_str.empty()){
            _params.gamma=parse_array<double>(gamma_str);
            for(size_t i=0; i<_params.gamma.size(); i++){
                _params.gamma[i]/=au2eV;
            }
        }
        else{
            throw std::string("gamma parameter is not specified!");
        }

        std::string s_str=ini.get("system").get("s");
        if(!s_str.empty()){
            _params.s=parse_array<double>(s_str);
        }
        else{
            throw std::string("s parameter is not specified!");
        }

        std::string Z_str=ini.get("system").get("Z");
        _params.Z=std::stod(Z_str);

        std::string Td_str=ini.get("system").get("Td");
        _params.Td=std::stod(Td_str)/au2fs;

        //parse output
        _params.outdir=ini.get("output").get("outdir");

        fs::path kgfile_path(_params.outdir);
        _params.kgfile_fname=(kgfile_path/=ini.get("output").get("kgfile")).c_str();

        fs::path rgfile_path(_params.outdir);
        _params.rgfile_fname=(rgfile_path/=ini.get("output").get("rgfile")).c_str();

        fs::path tfile_path(_params.outdir);
        _params.tfile_fname=(tfile_path/=ini.get("output").get("tfile")).c_str();

        fs::path densfile_path(_params.outdir);
        _params.densfile_fname=(densfile_path/=ini.get("output").get("densfile")).c_str();

        fs::path rhofile_path(_params.outdir);
        _params.rhofile_fname=(rhofile_path/=ini.get("output").get("rhofile")).c_str();

        fs::path afile_path(_params.outdir);
        _params.afile_fname=(afile_path/=ini.get("output").get("afile")).c_str();        

        fs::path pkfile_path(_params.outdir);
        _params.pkfile_fname=(pkfile_path/=ini.get("output").get("pkfile")).c_str();

        fs::path prfile_path(_params.outdir);
        _params.prfile_fname=(prfile_path/=ini.get("output").get("prfile")).c_str();        
    }
private:
    template<typename T>
    std::vector<T> parse_array(const std::string& str) const{
        std::vector<T> res;

        std::regex regex_number("[+-]?(\\d+([.]\\d*)?([eE][+-]?\\d+)?|[.]\\d+([eE][+-]?\\d+)?)");
        for(auto it=std::sregex_iterator(str.begin(),str.end(),regex_number);
                 it!=std::sregex_iterator(); it++){
            std::smatch m=*it;
            std::string valstr=m.str();

            std::istringstream ss(valstr);
            T val;
            ss>>val;

            res.push_back(val);
        }

        return res;
    }

    Parameters& _params;
};

#endif