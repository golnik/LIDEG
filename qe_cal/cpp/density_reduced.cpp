#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;
using namespace Eigen;
using namespace std;
using namespace std::chrono;

using Complex = std::complex<double>;
using RealVector = Vector3d;
using ComplexVector = vector<Complex>;
using RealMatrix = MatrixXd;

const double a = 2.46;

inline size_t idx_3d_to_1d(size_t i, size_t j, size_t k, size_t Nx, size_t Ny, size_t Nz)
{
    return i * Ny * Nz + j * Nz + k;
}

int extractNumber(const string &filename)
{
    size_t pos = filename.find_first_of("0123456789");
    return (pos == string::npos) ? -1 : stoi(filename.substr(pos));
}

void extract_kpoints_and_coefficients(const string &hdf5_directory, vector<double> &kpoints, vector<int> &g_vectors, vector<int> &g_size, vector<double> &coeffs_real_imag, int nbnd, size_t &total_kpoints)
{
    vector<string> coeff_files;
    for (const auto &entry : fs::directory_iterator(hdf5_directory))
    {
        if (entry.path().string().find("coeff_wfc") != string::npos)
        {
            coeff_files.push_back(entry.path().string());
        }
    }
    sort(coeff_files.begin(), coeff_files.end(), [](const string &a, const string &b)
         { return extractNumber(a) < extractNumber(b); });

    for (const auto &coeff_file : coeff_files)
    {
        ifstream infile(coeff_file);
        if (!infile)
        {
            cerr << "Error opening file: " << coeff_file << endl;
            continue;
        }

        RealVector k_vec;
        vector<int> g_vector(3);
        ComplexVector coeff_band(nbnd);

        string line;
        int current_band = -1;
        int g_coeff_size = 0;

        while (getline(infile, line))
        {
            if (line.find("xk:") != string::npos)
            {
                sscanf(line.c_str(), "xk: [%lf  %lf  %lf]", &k_vec[0], &k_vec[1], &k_vec[2]);
                kpoints.push_back(k_vec[0]);
                kpoints.push_back(k_vec[1]);
                kpoints.push_back(k_vec[2]);
            }
            else if (line.find("Total number of G vectors:") != string::npos)
            {
                sscanf(line.c_str(), "Total number of G vectors: %d", &g_coeff_size);
                g_size.push_back(g_coeff_size);
            }
            else if (line.find("band_index:") != string::npos)
            {
                sscanf(line.c_str(), "band_index: %d", &current_band);
                --current_band;
            }
            else if (line.find("[") != string::npos && current_band >= 0 && line.find("G vector") == string::npos)
            {
                double real_part, imag_part;
                sscanf(line.c_str(), "[%d, %d, %d] %lf %lf", &g_vector[0], &g_vector[1], &g_vector[2], &real_part, &imag_part);

                g_vectors.push_back(g_vector[0]);
                g_vectors.push_back(g_vector[1]);
                g_vectors.push_back(g_vector[2]);

                coeffs_real_imag.push_back(real_part);
                coeffs_real_imag.push_back(imag_part);
            }
        }
        ++total_kpoints;
    }
}

Complex compute_wavefunction_at_r(const double *k_vec, const int *g_vectors, vector<int> &g_size, const double *coeff_real_imag, const RealVector &r_vec)
{
    Complex wavefunction = 0.0;
    for (size_t j = 0; j < g_size.size(); j++)
    {
        for (size_t i = 0; i < g_size[i]; ++i)
        {
            RealVector g_new(
                (g_vectors[3 * i + 0] + g_vectors[3 * i + 1]) / sqrt(3),
                (g_vectors[3 * i + 0] - g_vectors[3 * i + 1]),
                g_vectors[3 * i + 2] * 0.245963);

            RealVector g_plus_k = Eigen::Map<const Vector3d>(k_vec) * (2 * M_PI / a) / 1.35177 + g_new;
            double phase = g_plus_k.dot(r_vec);

            Complex phase_factor = exp(Complex(0, 1) * phase);
            wavefunction += Complex(coeff_real_imag[2 * i], coeff_real_imag[2 * i + 1]) * phase_factor;
        }
    }

    return wavefunction;
}

MatrixXd compute_charge_density(const vector<double> &kpoints, const vector<int> &g_vectors, vector<int> &g_size, const vector<double> &coeffs_real_imag, const vector<RealVector> &real_space_grid, int nbnd, size_t total_kpoints)
{
    size_t grid_size = real_space_grid.size();
    MatrixXd rho_n_r(grid_size, nbnd);
    rho_n_r.setZero();

    std::cout << grid_size << std::endl;

#pragma omp parallel for collapse(2) // Enable parallelism with collapse for nested loops
    for (size_t k_idx = 0; k_idx < total_kpoints; ++k_idx)
    {
        for (int band_idx = 0; band_idx < nbnd; ++band_idx)
        {
            const double *k_vec = &kpoints[k_idx * 3];
            const int *g_vectors_band;
            const double *coeff_real_imag_band;
            if (k_idx == 0)
            {
                g_vectors_band = &g_vectors[k_idx * nbnd * 3];
                coeff_real_imag_band = &coeffs_real_imag[k_idx * nbnd * 2];
            }
            else
            {
                g_vectors_band = &g_vectors[k_idx * nbnd * 3 * g_size[k_idx - 1]];
                coeff_real_imag_band = &coeffs_real_imag[k_idx * nbnd * 2 * g_size[k_idx - 1]];
            }

#pragma omp parallel for
            for (size_t r_idx = 0; r_idx < grid_size; ++r_idx) //            for (size_t r_idx = 0; r_idx < grid_size; ++r_idx)
            {
                Complex psi_n_k_r = compute_wavefunction_at_r(k_vec, g_vectors_band, g_size, coeff_real_imag_band, real_space_grid[r_idx]);
#pragma omp atomic
                rho_n_r(r_idx, band_idx) += norm(psi_n_k_r) / total_kpoints;
            }
        }
    }

    // Print the progress after the parallel region
    for (size_t k_idx = 0; k_idx < total_kpoints; ++k_idx)
    {
        std::cout << k_idx << std::endl;
    }

    return rho_n_r;
}

void save_charge_density_to_file(const MatrixXd &rho_n_r, const vector<RealVector> &real_space_grid, const string &output_folder)
{
    if (!fs::exists(output_folder))
    {
        fs::create_directory(output_folder);
    }

#pragma omp parallel for
    for (int band_idx = 0; band_idx < rho_n_r.cols(); ++band_idx)
    {
        string filename = output_folder + "/rho_band_" + to_string(band_idx + 1) + ".dat";
        ofstream outfile(filename);
        if (!outfile)
        {
            cerr << "Error creating file: " << filename << endl;
            continue;
        }

        outfile << "x y z rho_n(r)\n";
        for (size_t r_idx = 0; r_idx < real_space_grid.size(); ++r_idx)
        {
            outfile << real_space_grid[r_idx][0] << " " << real_space_grid[r_idx][1] << " " << real_space_grid[r_idx][2] << " " << rho_n_r(r_idx, band_idx) << "\n";
        }
    }
    cout << "Charge density saved in " << output_folder << endl;
}

int main()
{
    auto start_total = high_resolution_clock::now();

    string hdf5_directory = "../coeff";
    string output_folder = "../wf";

    int nbnd = 8;
    int unit_cells_x = 1;
    int unit_cells_y = 1;
    int unit_cells_z = 1;
    int N0 = 8;

    vector<Vector3d> real_space_grid;
    Vector3d a1(unit_cells_x * (sqrt(3) / 2 * a), unit_cells_x * (1.0 / 2 * a), 0.0);
    Vector3d a2(unit_cells_y * (sqrt(3) / 2 * a), unit_cells_y * (-1.0 / 2 * a), 0.0);
    Vector3d a3(0.0, 0.0, unit_cells_z * 10.0);

    for (int i = 0; i < N0; ++i)
    {
        for (int j = 0; j < N0; ++j)
        {
            for (int k = 0; k < N0; ++k)
            {
                double frac_x = i / (double)(N0 - 1);
                double frac_y = j / (double)(N0 - 1);
                double frac_z = k / (double)(N0 - 1);

                Vector3d real_space_point = frac_x * a1 + frac_y * a2 + frac_z * a3;
                real_space_grid.push_back(real_space_point);
            }
        }
    }

    vector<double> kpoints;
    vector<int> g_vectors;
    vector<int> g_size;
    vector<double> coeffs_real_imag;
    size_t total_kpoints = 0;

    extract_kpoints_and_coefficients(hdf5_directory, kpoints, g_vectors, g_size, coeffs_real_imag, nbnd, total_kpoints);

    MatrixXd rho_n_r = compute_charge_density(kpoints, g_vectors, g_size, coeffs_real_imag, real_space_grid, nbnd, total_kpoints);

    save_charge_density_to_file(rho_n_r, real_space_grid, output_folder);

    auto end_total = high_resolution_clock::now();
    auto duration_total = duration_cast<milliseconds>(end_total - start_total);
    cout << "Total execution time: " << duration_total.count() << " milliseconds" << endl;

    return 0;
}
