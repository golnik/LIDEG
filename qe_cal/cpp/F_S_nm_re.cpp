#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;
using namespace Eigen;
using namespace std;
using namespace std::chrono; // For time management

// Type aliases for readability
using Complex = std::complex<double>;
using RealVector = Vector3d;
using ComplexVector = vector<Complex>;

// Lattice constant for graphene
const double a = 2.46;

// Helper function to extract the number from a file name
int extractNumber(const string &filename)
{
    size_t pos = filename.find_first_of("0123456789");
    return (pos == string::npos) ? -1 : stoi(filename.substr(pos));
}

// Step 1: Extract k-point information and coefficients from .dat files
#include <cmath> // For sqrt

void extract_kpoints_and_coefficients(const string &hdf5_directory, vector<RealVector> &kpoints, vector<vector<vector<vector<int>>>> &g_vectors_list, vector<vector<ComplexVector>> &coeffs_real_imag, int nbnd)
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
        vector<vector<vector<int>>> g_vectors(nbnd); // Initialize per band, store vectors of G-vectors (3D integer vectors)
        vector<ComplexVector> coeff_real_imag(nbnd);

        string line;
        int current_band = -1;
        while (getline(infile, line))
        {
            if (line.find("xk:") != string::npos)
            {
                sscanf(line.c_str(), "xk: [%lf %lf %lf]", &k_vec[0], &k_vec[1], &k_vec[2]);
                // std::cout << k_vec << std::endl;
                kpoints.push_back(k_vec);
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

                g_vectors[current_band].push_back(g_vector);
                coeff_real_imag[current_band].emplace_back(real_part, imag_part);
            }
        }

        // Normalize coefficients for each band
        for (int band = 0; band < nbnd; ++band)
        {
            // Compute the sum of |coeff_real_imag|^2 for the current band
            double norm_factor = 0.0;
            for (const auto &coeff : coeff_real_imag[band])
            {
                double real_part = coeff.real();
                double imag_part = coeff.imag();
                norm_factor += real_part * real_part + imag_part * imag_part;
            }

            if (norm_factor > 0.0)
            {
                norm_factor = sqrt(norm_factor); // Take the square root to get the normalization factor
                // std::cout << norm_factor << std::endl;
                // Normalize each coefficient
                for (auto &coeff : coeff_real_imag[band])
                {
                    coeff /= norm_factor; // Normalize the real and imaginary parts
                }
            }
        }

        g_vectors_list.push_back(g_vectors);
        coeffs_real_imag.push_back(coeff_real_imag);
    }
}

void rearrange_kpoints(const vector<RealVector> &kpoints, const vector<vector<vector<vector<int>>>> &g_vectors_list, const vector<vector<ComplexVector>> &coeffs_real_imag, int nbnd, int m_min, int n_min)
{
    RealVector S(1, 0, 0); // Vector S = [1, 0, 0]

    // Find min/max kx and ky
    double min_kx = std::numeric_limits<double>::max();
    double max_kx = std::numeric_limits<double>::lowest();
    double min_ky = std::numeric_limits<double>::max();
    double max_ky = std::numeric_limits<double>::lowest();

    RealVector xL, xR; // xL for min kx, xR for max kx
    RealVector yL, yR; // yL for min ky, yR for max ky

    // Find min and max kx, ky (sequential part)
    for (const auto &k_vec : kpoints)
    {
        if (k_vec[0] < min_kx)
        {
            min_kx = k_vec[0];
            xL = k_vec;
        }
        if (k_vec[0] > max_kx)
        {
            max_kx = k_vec[0];
            xR = k_vec;
        }
        if (k_vec[1] < min_ky)
        {
            min_ky = k_vec[1];
            yL = k_vec;
        }
        if (k_vec[1] > max_ky)
        {
            max_ky = k_vec[1];
            yR = k_vec;
        }
    }

    double b0 = xR[0] - xL[0];
    int total_kpoints = kpoints.size();
    int x = static_cast<int>(sqrt(total_kpoints));

    std::cout << xR[0] - xL[0] << std::endl;
    std::cout << yR[1] - yL[1] << std::endl;
    std::cout << x << std::endl;

    if (x * x != total_kpoints)
    {
        cerr << "Error: The number of k-points is not a perfect square." << endl;
        return;
    }

    // Precompute reciprocal space basis vectors
    RealVector b1 = RealVector(b0 / 2, b0 * sqrt(3) / 2, 0.0) / (x - 1);
    RealVector b2 = RealVector(b0 / 2, -b0 * sqrt(3) / 2, 0.0) / (x - 1);

    std::cout << b1 << std::endl;
    std::cout << b2 << std::endl;

    string output_dir = "../transition_re";
    if (!fs::exists(output_dir))
    {
        fs::create_directory(output_dir);
    }

    // Use all available threads (performance + efficiency cores)
    omp_set_num_threads(16); // Adjust based on your preference or hyper-threading efficiency

// Parallelize the outer loop over kpoints using dynamic scheduling
#pragma omp parallel for schedule(dynamic)
    for (int old_idx = 0; old_idx < total_kpoints; ++old_idx)
    {
        RealVector k_vec = kpoints[old_idx] - xL; // Subtract origin O

        // Solve for N1 and N2
        Matrix2d basis_matrix;
        basis_matrix << b1[0], b2[0], b1[1], b2[1];
        Vector2d N_vals = basis_matrix.colPivHouseholderQr().solve(k_vec.head<2>());

        int N1 = static_cast<int>(round(N_vals[0]));
        int N2 = static_cast<int>(round(N_vals[1]));

        // Compute the unique new index based on N1 and N2
        int new_index = N1 * x + N2;

        // Create output file (use thread-safe filename generation)
        string filename = output_dir + "/" + to_string(new_index) + "_rearrange.dat";
        ofstream outfile(filename, ios::out); // Open the output file

        // Write k-point and other data to file
        outfile << "xk: [" << kpoints[old_idx].transpose() << "]" << std::endl;
        outfile << "S : [" << S[0] << '\t' << S[1] << '\t' << S[2] << "]" << std::endl;

        for (int band_m = m_min - 1; band_m < nbnd; ++band_m)
        {
            for (int band_n = max(n_min - 1, band_m); band_n < nbnd; ++band_n)
            {
                Complex sum = 0.0;

                const auto &g_vectors_m = g_vectors_list[old_idx][band_m];
                const auto &g_vectors_n = g_vectors_list[old_idx][band_n];
                const auto &coeff_m = coeffs_real_imag[old_idx][band_m];
                const auto &coeff_n = coeffs_real_imag[old_idx][band_n];

                // Compute G1, G2 and sum for bands
                for (size_t i = 0; i < g_vectors_m.size(); ++i)
                {
                    for (size_t j = 0; j < g_vectors_n.size(); ++j)
                    {
                        RealVector G1(g_vectors_m[i][0], g_vectors_m[i][1], g_vectors_m[i][2]);
                        RealVector G2(g_vectors_n[j][0], g_vectors_n[j][1], g_vectors_n[j][2]);

                        if ((S + G1 - G2).norm() < 1e-6) // Floating-point comparison
                        {
                            sum += conj(coeff_m[i]) * coeff_n[j];
                        }
                    }
                }

                // Write results for this (m, n) pair
                outfile << "Sum for band m = " << band_m + 1 << ", band n = " << band_n + 1 << ": ("
                        << real(sum) << "," << imag(sum) << ")\n";
            }
        }

        outfile.close(); // Close the output file

// Output status (optional, but can be disabled for performance)
#pragma omp critical
        {
            cout << "Results for k-point " << old_idx + 1 << " saved in " << filename << endl;
        }
    }
}

int main()
{
    // Define input parameters
    string hdf5_directory = "../coeff"; // Path to directory containing coeff_wfcXXX.dat files
    int nbnd = 8;                       // Number of bands
    int m_min = 4;                      // Minimum value of m (condition m >= m_min)
    int n_min = 4;                      // Minimum value of n (condition n >= m and n >= n_min)

    // Extract kpoints and coefficients (not parallelized, since it's I/O-bound)
    vector<RealVector> kpoints;
    vector<vector<vector<vector<int>>>> g_vectors_list;
    vector<vector<ComplexVector>> coeffs_real_imag;

    extract_kpoints_and_coefficients(hdf5_directory, kpoints, g_vectors_list, coeffs_real_imag, nbnd);

    // Rearrange kpoints and perform computations in parallel
    rearrange_kpoints(kpoints, g_vectors_list, coeffs_real_imag, nbnd, m_min, n_min);

    return 0;
}