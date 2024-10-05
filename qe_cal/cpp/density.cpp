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
using namespace std::chrono; // For time management

// Type aliases for readability
using Complex = std::complex<double>;
using RealVector = Vector3d;
using ComplexVector = vector<Complex>;
using RealMatrix = MatrixXd;

double a = 2.46; // Lattice constant for graphene

// Helper function to extract the number from a file name
int extractNumber(const string &filename)
{
    size_t pos = filename.find_first_of("0123456789");
    return (pos == string::npos) ? -1 : stoi(filename.substr(pos));
}

// Step 1: Extract k-point information and coefficients from .dat files
void extract_kpoints_and_coefficients(const string &hdf5_directory, vector<RealVector> &kpoints, vector<vector<vector<vector<int>>>> &g_vectors_list, vector<vector<ComplexVector>> &coeffs_real_imag, int nbnd)
{
    // List the files and sort by number
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
        vector<vector<vector<int>>> g_vectors(nbnd); // Initialize per band, store vectors of G-vectors (3D integer vectors)
        vector<ComplexVector> coeff_real_imag(nbnd);

        string line;
        int current_band = -1;
        while (getline(infile, line))
        {
            if (line.find("xk:") != string::npos)
            {
                sscanf(line.c_str(), "xk: [%lf %lf %lf]", &k_vec[0], &k_vec[1], &k_vec[2]);
                kpoints.push_back(k_vec); //                kpoints.push_back(k_vec);
            }
            else if (line.find("band_index:") != string::npos)
            {
                sscanf(line.c_str(), "band_index: %d", &current_band);
                --current_band; // Convert to 0-indexed
            }
            else if (line.find("[") != string::npos && current_band >= 0 && line.find("G vector") == string::npos)
            {
                vector<int> g_vector(3);
                double real_part, imag_part;
                sscanf(line.c_str(), "[%d, %d, %d] %lf %lf", &g_vector[0], &g_vector[1], &g_vector[2], &real_part, &imag_part);

                g_vectors[current_band].push_back(g_vector);                      // Append G-vector to the current band's list
                coeff_real_imag[current_band].emplace_back(real_part, imag_part); // Append the complex coefficient
            }
        }

        g_vectors_list.push_back(g_vectors);
        coeffs_real_imag.push_back(coeff_real_imag);
    }
}

// Step 2: Compute the wavefunction for each k-point and band at a given real-space point
Complex compute_wavefunction_at_r(const RealVector &k_vec, const vector<vector<int>> &g_vectors, const ComplexVector &coeff_real_imag, const RealVector &r_vec)
{
    Complex wavefunction = 0.0;
    // std::cout << g_vectors.size() << std::endl;
    // psi_{n,k}(r) = sum c_i * exp(i * (G_i + k) · r)
    for (size_t i = 0; i < g_vectors.size(); ++i) //    for (size_t i = 0; i < g_vectors.size(); ++i)
    {

        const auto &g = g_vectors[i];
        RealVector g_new(
            (g[0] + g[1]) / sqrt(3), // First component
            (g[0] - g[1]),           // Second component
            g[2] * 0.245963          // Third component
        );
        // std::cout << k_vec << std::endl;
        // Now scale by (2 * M_PI / a)
        g_new *= (2 * M_PI / a); // Proper way to multiply the whole vector by a scalar

        // Scale k_vec by 1.35177 and add g_new
        RealVector g_plus_k = k_vec * 0 + g_new;//        RealVector g_plus_k = k_vec * (2 * M_PI / a) / 1.35177 + g_new;

        double phase = g_plus_k.dot(r_vec);
        Complex phase_factor = exp(Complex(0, 1) * phase);
        wavefunction += coeff_real_imag[i] * phase_factor;
    }

    return wavefunction;
}

// Step 3: Compute rho_n(r) by summing |psi_{n,k}(r)|^2 over all k-points
MatrixXd compute_charge_density(const vector<RealVector> &kpoints, const vector<vector<vector<vector<int>>>> &g_vectors_list, const vector<vector<ComplexVector>> &coeffs_real_imag, const vector<RealVector> &real_space_grid, int nbnd)
{
    size_t total_kpoints = kpoints.size();
    size_t grid_size = real_space_grid.size();
    MatrixXd rho_n_r(grid_size, nbnd);
    rho_n_r.setZero();

    // Parallelize the outer loops using OpenMP
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (size_t k_idx = 0; k_idx < total_kpoints; ++k_idx)
    {
        for (int band_idx = 0; band_idx < nbnd; ++band_idx)
        {
            const RealVector &k_vec = kpoints[k_idx];
            const auto &g_vectors = g_vectors_list[k_idx][band_idx];
            const auto &coeff_real_imag = coeffs_real_imag[k_idx][band_idx];

            // Step 1: Compute normalization factor
            double norm_factor = 0.0;

            // Compute the total sum of |psi_n_k_r|^2 for normalization
#pragma omp parallel for reduction(+:norm_factor) schedule(static)
            for (size_t r_idx = 0; r_idx < grid_size; ++r_idx)
            {
                Complex psi_n_k_r = compute_wavefunction_at_r(k_vec, g_vectors, coeff_real_imag, real_space_grid[r_idx]);
                norm_factor += norm(psi_n_k_r); // Accumulate |psi|^2
            }

            // Normalize the factor by the total number of points
            norm_factor = sqrt(norm_factor); // Take the square root to get the norm

            // Step 2: Recompute wavefunction contributions and update rho_n_r
#pragma omp parallel for schedule(static)
            for (size_t r_idx = 0; r_idx < grid_size; ++r_idx)
            {
                Complex psi_n_k_r = compute_wavefunction_at_r(k_vec, g_vectors, coeff_real_imag, real_space_grid[r_idx]);

                // Normalize psi_n_k_r by norm_factor
                psi_n_k_r /= norm_factor;

#pragma omp atomic
                rho_n_r(r_idx, band_idx) += norm(psi_n_k_r) / total_kpoints; // norm gives |psi|^2
            }
        }
    }

    return rho_n_r;
}


// Step 4: Save rho_n(r) to file
void save_charge_density_to_file(const MatrixXd &rho_n_r, const vector<RealVector> &real_space_grid, const string &output_folder)
{
    if (!fs::exists(output_folder))
    {
        fs::create_directory(output_folder);
    }

#pragma omp parallel for // Parallelize the saving of each band
    for (int band_idx = 0; band_idx < rho_n_r.cols(); ++band_idx)
    {
        string filename = output_folder + "/rho_band_" + to_string(band_idx + 1) + ".dat";
        ofstream outfile(filename);
        if (!outfile)
        {
            cerr << "Error creating file: " << filename << endl;
            continue;
        }

        outfile << "x y z rho_n(r)\n"; // Now include z coordinate
        for (size_t r_idx = 0; r_idx < real_space_grid.size(); ++r_idx)
        {
            outfile << real_space_grid[r_idx][0] << " " << real_space_grid[r_idx][1] << " " << real_space_grid[r_idx][2] << " " << rho_n_r(r_idx, band_idx) << "\n";
        }
    }
    cout << "Charge density saved in " << output_folder << endl;
}

int main()
{
    // Start timer for the entire process
    auto start_total = high_resolution_clock::now();

    string hdf5_directory = "../coeff"; // Path to directory containing coeff_wfcXXX.dat files
    string output_folder = "../wf";

    int nbnd = 8;

    int unit_cells_x = 1; // Number of unit cells along x
    int unit_cells_y = 1; // Number of unit cells along y
    int unit_cells_z = 1; // Number of unit cells along z (for 2D graphene, we often keep z = 1)

    // Keep the number of points per unit cell constant
    int Nx = 16, Ny = 16, Nz = 16; // These represent the number of grid points per unit cell

    vector<Vector3d> real_space_grid;

    // Define the graphene unit cell lattice vectors and scale them by the number of unit cells
    Vector3d a1(unit_cells_x * (sqrt(3) / 2 * a), unit_cells_x * (1.0 / 2 * a), 0.0);  // First lattice vector scaled by unit_cells_x
    Vector3d a2(unit_cells_y * (sqrt(3) / 2 * a), unit_cells_y * (-1.0 / 2 * a), 0.0); // Second lattice vector scaled by unit_cells_y
    Vector3d a3(0.0, 0.0, unit_cells_z * 10.0);                                        // New vector along the z-axis, scaled by unit_cells_z

    // Start timer for grid generation
    auto start_grid = high_resolution_clock::now();

    // Loop to generate real-space grid points, now including the z direction
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int k = 0; k < Nz; ++k) // Adding z-layer grid loop
            {
                // Adjust the fractional coordinates to include boundary points
                double frac_x = i / (double)(Nx - 1); // Ensure that frac_x ranges from 0 to 1 inclusive
                double frac_y = j / (double)(Ny - 1); // Ensure that frac_y ranges from 0 to 1 inclusive
                double frac_z = k / (double)(Nz - 1); // Ensure that frac_z ranges from 0 to 1 inclusive

                // Convert fractional coordinates to real-space coordinates using scaled lattice vectors
                Vector3d real_space_point = frac_x * a1 + frac_y * a2 + frac_z * a3;

                // Store the real-space point
                real_space_grid.push_back(real_space_point);
            }
        }
    }

    auto end_grid = high_resolution_clock::now();
    auto duration_grid = duration_cast<milliseconds>(end_grid - start_grid);
    // cout << "Grid generation time: " << duration_grid.count() << " milliseconds" << endl;

    vector<RealVector> kpoints;
    vector<vector<vector<vector<int>>>> g_vectors_list;
    vector<vector<ComplexVector>> coeffs_real_imag;

    // Start timer for k-point extraction
    auto start_kpoints = high_resolution_clock::now();

    // Extract kpoints and coefficients
    extract_kpoints_and_coefficients(hdf5_directory, kpoints, g_vectors_list, coeffs_real_imag, nbnd);

    auto end_kpoints = high_resolution_clock::now();
    auto duration_kpoints = duration_cast<milliseconds>(end_kpoints - start_kpoints);
    // cout << "K-point extraction time: " << duration_kpoints.count() << " milliseconds" << endl;

    // Start timer for charge density computation
    auto start_density = high_resolution_clock::now();

    // Compute charge density
    MatrixXd rho_n_r = compute_charge_density(kpoints, g_vectors_list, coeffs_real_imag, real_space_grid, nbnd);

    auto end_density = high_resolution_clock::now();
    auto duration_density = duration_cast<milliseconds>(end_density - start_density);
    // cout << "Charge density computation time: " << duration_density.count() << " milliseconds" << endl;

    // Start timer for file saving
    auto start_save = high_resolution_clock::now();

    // Save charge density
    save_charge_density_to_file(rho_n_r, real_space_grid, output_folder);

    auto end_save = high_resolution_clock::now();
    auto duration_save = duration_cast<milliseconds>(end_save - start_save);
    // cout << "File saving time: " << duration_save.count() << " milliseconds" << endl;

    // Total time
    auto end_total = high_resolution_clock::now();
    auto duration_total = duration_cast<milliseconds>(end_total - start_total);
    // cout << "Total execution time: " << duration_total.count() << " milliseconds" << endl;

    return 0;
}
