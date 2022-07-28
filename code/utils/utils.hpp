#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <complex>

#include "integrator.hpp"

#define au2nm 0.052917829614246
#define au2A  0.52917829614246
#define au2eV 27.211396641308
#define au2fs 0.02418884254

typedef std::complex<double> complex_t;
const complex_t I=complex_t(0.,1.);

#include <boost/algorithm/string.hpp>

#include <algorithm> 
#include <cctype>
#include <locale>

enum kgrid_types {regular=0, quad};

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

void read_column_from_file(const std::string& fname, const size_t& icol,
    std::vector<double>& data){
    
    std::ifstream inpstream(fname);

    std::string line;
    while(std::getline(inpstream,line)){
        if(line.front()!='#'){
            std::vector<std::string> split_str;
            std::istringstream iss(line);
            for(std::string s; iss>>s;) 
                split_str.push_back(s);

            double value=std::stod(split_str[icol]);
            data.push_back(value);
        }
    }

    return;
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::vector<double> create_grid(const double& min, const double& max, const size_t& N, const int& type=0){
    std::vector<double> grid(N);

    if(type==regular){
        double d=(max-min)/(N-1);
        double var=min;
        for(size_t i=0; i<N; i++){
            grid[i]=var;
            var+=d;
        }
    }
    else if(type==quad){
        grid=create_gauss_legendre_grid(min,max,N);
    }

    return grid;
}

template<typename matrix_t>
matrix_t read_2D_from_file(const std::string& fname, const size_t& icol, const size_t& N, const size_t& M){
    matrix_t res(N,M);

    std::vector<double> data;
    read_column_from_file(fname,icol,data);

    int indx=0;
    for(size_t i=0; i<N; i++){
        for(size_t j=0; j<M; j++){
            //std::cout<<indx<<" "<<data[indx]<<std::endl;
            res(i,j)=data[indx];
            indx++;
        }
    }

    return res;
}

#endif