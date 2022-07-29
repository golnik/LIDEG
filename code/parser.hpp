#ifndef PARSER_HPP
#define PARSER_HPP

#include "mini/ini.h"
#include "utils/utils.hpp"

#include <string>
#include <fstream>

#include <filesystem>
namespace fs=std::filesystem;

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

    int kgrid_type;

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

    double E0;
    std::string field_fname;

    double a;
    double e2p;
    double gamma;
    double s;
    double Z;

    std::string outdir;
    std::string gfile_fname;
    std::string tfile_fname;
    std::string densfile_fname;
    std::string rhofile_fname;
    std::string afile_fname;

    void print(std::ostream& out) const{
        out<<tmin<<" "<<tmax<<" "<<Nt<<std::endl;
        out<<err_abs<<" "<<err_rel<<std::endl;

        //out<<"kx grid:"<<std::endl;
        //out<<kx_min<<" "<<kx_max<<" "<<Nkx<<std::endl;
        //out<<"ky grid:"<<std::endl;
        //out<<ky_min<<" "<<ky_max<<" "<<Nky<<std::endl;

        out<<field_fname<<std::endl;

        out<<"Parameters of the model:"<<std::endl;
        out<<"a: "<<a<<std::endl;
        out<<"e2p: "<<e2p<<std::endl;
        out<<"gamma: "<<gamma<<std::endl;
        out<<"s: "<<s<<std::endl;
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
        if(kgrid_type_str=="regular"){
            _params.kgrid_type=regular;
        }
        else if(kgrid_type_str=="quad"){
            _params.kgrid_type=quad;
        }
        else{
            throw std::string("Wrong kgrid_type!");
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

        //parse field file
        _params.E0=std::stod(ini.get("field").get("E0"))/au2Vnm;
        _params.field_fname=ini.get("field").get("fname");

        //parse model parameters
        std::string a_str=ini.get("model").get("a");
        _params.a=std::stod(a_str)/au2A;

        std::string e2p_str=ini.get("model").get("e2p");
        _params.e2p=std::stod(e2p_str)/au2eV;

        std::string gamma_str=ini.get("model").get("gamma");
        _params.gamma=std::stod(gamma_str)/au2eV;

        std::string s_str=ini.get("model").get("s");
        _params.s=std::stod(s_str);

        std::string Z_str=ini.get("model").get("Z");
        _params.Z=std::stod(Z_str);

        //parse output
        _params.outdir=ini.get("output").get("outdir");

        fs::path gfile_path(_params.outdir);
        _params.gfile_fname=(gfile_path/=ini.get("output").get("gfile")).c_str();

        fs::path tfile_path(_params.outdir);
        _params.tfile_fname=(tfile_path/=ini.get("output").get("tfile")).c_str();

        fs::path densfile_path(_params.outdir);
        _params.densfile_fname=(densfile_path/=ini.get("output").get("densfile")).c_str();

        fs::path rhofile_path(_params.outdir);
        _params.rhofile_fname=(rhofile_path/=ini.get("output").get("rhofile")).c_str();

        fs::path afile_path(_params.outdir);
        _params.afile_fname=(afile_path/=ini.get("output").get("afile")).c_str();        
    }
private:
    Parameters& _params;
};

#endif