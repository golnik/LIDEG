#ifndef PARSER_HPP
#define PARSER_HPP

#include "mini/ini.h"
#include "utils.hpp"

#include <string>
#include <fstream>

struct Parameters{
    double tmin;
    double tmax;
    size_t Nt;

    double ddt;
    double err_abs;
    double err_rel;

    double kx_min;
    double kx_max;
    size_t Nkx;

    double ky_min;
    double ky_max;
    size_t Nky;

    size_t Nclx;
    size_t Ncly;

    double xmin;
    double xmax;
    size_t Nx;

    double ymin;
    double ymax;
    size_t Ny;

    double zmin;
    double zmax;
    size_t Nz;        

    std::string field_fname;

    double a;
    double e2p;
    double gamma;
    double s;
    double Z;

    std::string tfile_fname;
    std::string rhofile_fname;
    
    std::string densfile_fname;

    void print(std::ostream& out) const{
        out<<tmin<<" "<<tmax<<" "<<Nt<<std::endl;
        out<<err_abs<<" "<<err_rel<<std::endl;

        out<<"kx grid:"<<std::endl;
        out<<kx_min<<" "<<kx_max<<" "<<Nkx<<std::endl;
        out<<"ky grid:"<<std::endl;
        out<<ky_min<<" "<<ky_max<<" "<<Nky<<std::endl;

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
        std::string kx_min_str=ini.get("kgrid").get("kx_min");
        _params.kx_min=std::stod(kx_min_str)*au2nm;
        
        std::string kx_max_str=ini.get("kgrid").get("kx_max");
        _params.kx_max=std::stod(kx_max_str)*au2nm;

        std::string Nkx_str=ini.get("kgrid").get("Nkx");
        _params.Nkx=std::stoi(Nkx_str);  

        std::string ky_min_str=ini.get("kgrid").get("ky_min");
        _params.ky_min=std::stod(ky_min_str)*au2nm;

        std::string ky_max_str=ini.get("kgrid").get("ky_max");
        _params.ky_max=std::stod(ky_max_str)*au2nm;

        std::string Nky_str=ini.get("kgrid").get("Nky");
        _params.Nky=std::stoi(Nky_str);

        //parse rgrid section
        std::string Nclx_str=ini.get("rgrid").get("Nclx");
        _params.Nclx=std::stoi(Nclx_str);
        std::string Ncly_str=ini.get("rgrid").get("Ncly");
        _params.Ncly=std::stoi(Ncly_str);

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
        _params.tfile_fname=ini.get("output").get("tfile");
        _params.rhofile_fname=ini.get("output").get("rhofile");

        _params.densfile_fname=ini.get("output").get("densfile");
    }
private:
    Parameters& _params;
};

#endif