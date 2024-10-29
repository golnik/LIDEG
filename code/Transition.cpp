#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "mini/ini.h"

#include "external_field.hpp"
#include "graphenemodel.hpp"
#include "model1/graphene.hpp"
#include "parser.hpp"
#include "utils/grid.hpp"
#include "utils/utils.hpp"
// #include "model1/WFs.hpp"
#include "Nlayer/nlayer.hpp"
#include "model2/graphene2.hpp"

#include "WFs.hpp"

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;
typedef vector<complex_t> vector_t;
typedef matrix<double> matrix_t;

#include <boost/format.hpp>

int main(int argc, char **argv)
{
    try
    {
        std::string fname = argv[1];

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        // params.print(std::cout);

        // read fields from file
        std::vector<double> tgrid_fit;
        std::vector<double> Adata_x_fit;
        std::vector<double> Adata_y_fit;
        std::vector<double> Edata_x_fit;
        std::vector<double> Edata_y_fit;

        read_column_from_file(params.field_fname, 0, tgrid_fit);
        read_column_from_file(params.field_fname, 1, Adata_x_fit); /// should be x
        read_column_from_file(params.field_fname, 2, Adata_y_fit); /// should be y
        read_column_from_file(params.field_fname, 3, Edata_x_fit); /// should be x
        read_column_from_file(params.field_fname, 4, Edata_y_fit); /// should be y

        double t0_fit = tgrid_fit[0] / au2fs;
        double dt_fit = (tgrid_fit[1] - tgrid_fit[0]) / au2fs;

        ExternalField *E0 = nullptr;

        ExternalField *Afield_x = new ExternalFieldFromData(Adata_x_fit, t0_fit, dt_fit, params.E0);
        ExternalField *Afield_y = new ExternalFieldFromData(Adata_y_fit, t0_fit, dt_fit, params.E0);
        ExternalField *Efield_x = new ExternalFieldFromData(Edata_x_fit, t0_fit, dt_fit, params.E0);
        ExternalField *Efield_y = new ExternalFieldFromData(Edata_y_fit, t0_fit, dt_fit, params.E0);

        auto tgrid = create_grid(params.tmin, params.tmax, params.Nt);

        for (size_t tstep = 0; tstep < 1; tstep++)
        {
            double time = tgrid[tstep];

            // get vector potential at the given time step
            double Ax = (*Afield_x)(time);
            double Ay = (*Afield_y)(time);

            // std::cout << "Ax: " << Ax << " Ay: " << Ay << std::endl;

            /// prepare Black spots
            typedef std::vector<std::vector<double>> BZ_t;
            BZ_t BZ0{{0, 0}};
            BZ_t BZ1{{1, 0}, {1, 1}, {0, 1}, {-1, 0}, {-1, -1}, {0, -1}};
            BZ_t BZ2{{2, 1}, {1, 2}, {-1, 1}, {-2, -1}, {-1, -2}, {1, -1}};
            BZ_t BZ3{{2, 0}, {2, 2}, {0, 2}, {-2, 0}, {-2, -2}, {0, -2}};

            BZ_t BZt{{1, 0}};

            std::vector<BZ_t> zones{BZt};

            size_t nzones = zones.size();
            size_t nspots = 0;
            for (size_t izone = 0; izone < nzones; izone++)
            {
                for (auto spot : zones[izone])
                {
                    nspots++;
                }
            }

            ///////////////////

            // parameters

            //////////////////

            size_t Nx_qe = 20;
            size_t Ny_qe = 20;
            size_t Nz_qe = 80;

            size_t Nst_qe = 8;

            // create rgrid
            Grid2D *xygrid;
            if (params.rgrid_type == rgrid_types::rectan)
            {
                Grid1D *xgrid = new RegularGrid1D(params.xmin, params.xmax, Nx_qe);
                Grid1D *ygrid = new RegularGrid1D(params.ymin, params.ymax, Ny_qe);

                xygrid = new RegularGrid2D(xgrid, ygrid);
            }
            else if (params.rgrid_type == rgrid_types::ucell)
            {
                double Ox = -(1. / sqrt(3.)) * params.a;
                double Oy = 0.;
                double a1x = params.a / 2. * sqrt(3.);
                double a1y = params.a / 2.;
                double a2x = a1x;
                double a2y = -a1y;

                xygrid = new UCellGrid2D(Ox, Oy, a1x, a1y, Nx_qe, a2x, a2y, Ny_qe);
            }

            // std::cout << params.a << std::endl;

            Grid1D *zgrid = new RegularGrid1D(0, 18.895, Nz_qe);

            Integrator2D *integrator_xy = new Integrator2D(xygrid);
            Integrator1D *integrator_z = new Integrator1D(zgrid);

            MultiIndex indx_xyz({Nx_qe, Ny_qe, Nz_qe});
            size_t N_xyz = indx_xyz.size();

            /////////////////////////////////////

            MultiIndex indx_xy({Nx_qe, Ny_qe});
            size_t N_xy = indx_xy.size();

            // std::vector<matrix_t> dens_data(Nst_qe);

            //////////

            // The 1st col is time (indx 0)

            // For each Bragg spot start from 2nd col (indx 1)

            // indx [3n + 1] is intra
            // indx [3n + 2] is inter
            // indx [3n + 3] is total
            // where n is nature number
            size_t Nk = params.Nkx * params.Nky;
#pragma omp parallel for num_threads(90) schedule(dynamic)
            for (size_t ik = 1; ik < Nk + 1; ik++)
            {
                /* code */

                // Construct the output file name based on the k-point
                std::string filename = "/xdisk/ngolubev/mingruiyuan/QE_diffr/transition_re/" + std::to_string(ik) + "_rearrange.dat";
                std::ofstream outfile(filename, std::ios::out); // Open the output file
                if (!outfile.is_open())
                {
                    throw std::runtime_error("Unable to open output file: " + filename);
                }

                int t_indx = 0;
                size_t ispot = 0;
                for (size_t izone = 0; izone < nzones; izone++)
                {
                    for (auto spot : zones[izone])
                    {
                        auto m = spot[0];
                        auto n = spot[1];

                        for (size_t mst = 4; mst < Nst_qe + 1; mst++)
                        {
                            for (size_t fst = 4; fst < Nst_qe + 1; fst++)
                            {
                                /////////////////////////////////

                                // generate array for wf

                                /////////////////////////////////

                                /////////////////////////////////
                                // Read wavefunction data for mst
                                /////////////////////////////////

                                std::string file_name_mst = "/xdisk/ngolubev/mingruiyuan/QE_diffr/wfc/wfc_" + std::to_string(mst) + "_" + std::to_string(ik) + ".dat";
                                std::ifstream file_mst(file_name_mst);
                                if (!file_mst.is_open())
                                {
                                    throw std::runtime_error("File not found: " + file_name_mst);
                                }

                                std::vector<complex_t> psi_mst_k;
                                double real_part, imag_part;
                                while (file_mst >> real_part >> imag_part)
                                {
                                    psi_mst_k.emplace_back(real_part, imag_part);
                                }
                                file_mst.close();

                                // Check if the data size matches the expected dimensions
                                size_t expected_size = static_cast<size_t>(Nx_qe) * Ny_qe * Nz_qe;
                                if (psi_mst_k.size() != expected_size)
                                {
                                    std::cout << file_name_mst << std::endl;
                                    std::cout << psi_mst_k.size() << std::endl;

                                    throw std::runtime_error("Data size for mst does not match the expected dimensions.");
                                }

                                /////////////////////////////////
                                // Read wavefunction data for fst
                                /////////////////////////////////
                                std::string file_name_fst = "/xdisk/ngolubev/mingruiyuan/QE_diffr/wfc/wfc_" + std::to_string(fst) + "_" + std::to_string(ik) + ".dat";
                                std::ifstream file_fst(file_name_fst);
                                if (!file_fst.is_open())
                                {
                                    throw std::runtime_error("File not found: " + file_name_fst);
                                }

                                std::vector<complex_t> psi_fst_k;
                                while (file_fst >> real_part >> imag_part)
                                {
                                    psi_fst_k.emplace_back(real_part, imag_part);
                                }
                                file_fst.close();

                                // Check if the data size matches the expected dimensions
                                if (psi_fst_k.size() != expected_size)
                                {
                                    throw std::runtime_error("Data size for fst does not match the expected dimensions.");
                                }

                                ///////////////////////////////

                                // z int

                                ///////////////////////////////

                                std::vector<complex_t> Q_xy_f(N_xy);
                                std::vector<complex_t> Q_xy_m(N_xy);

                                /////////////////////

                                // normalization

                                /////////////////////

                                for (size_t ix = 0; ix < Nx_qe; ix++)
                                {
                                    for (size_t iy = 0; iy < Ny_qe; iy++)
                                    {
                                        std::complex<double> result = 0.0;
                                        integrator_z->trapz(
                                            [ix, iy, &psi_fst_k, Nx_qe, Ny_qe](const size_t &iz)
                                            {
                                                size_t index = iz * (Ny_qe * Nx_qe) + iy * Nx_qe + ix;
                                                return std::conj(psi_fst_k[index]) * psi_fst_k[index];
                                            },
                                            result);

                                        // Do something with the result, e.g., store it in a 2D array
                                        size_t indx_ixiy = indx_xy({iy, ix});
                                        Q_xy_f[indx_ixiy] = result;
                                    }
                                }

#pragma omp parallel for collapse(2) schedule(static)
                                for (size_t ix = 0; ix < Nx_qe; ix++)
                                {
                                    for (size_t iy = 0; iy < Ny_qe; iy++)
                                    {
                                        std::complex<double> result = 0.0;
                                        integrator_z->trapz(
                                            [ix, iy, &psi_mst_k, Nx_qe, Ny_qe](const size_t &iz)
                                            {
                                                size_t index = iz * (Ny_qe * Nx_qe) + iy * Nx_qe + ix;
                                                return std::conj(psi_mst_k[index]) * psi_mst_k[index];
                                            },
                                            result);

                                        // Do something with the result, e.g., store it in a 2D array
                                        size_t indx_ixiy = indx_xy({iy, ix});
                                        Q_xy_m[indx_ixiy] = result;
                                    }
                                }

                                auto density = [m, n, &xygrid, &indx_xy, &params](const size_t &iy, const size_t &ix, const std::vector<complex_t> &Q_xy)
                                {
                                    double x = (*xygrid)(ix, iy)[0];
                                    double y = (*xygrid)(ix, iy)[1];

                                    // size_t indx_ixiy = indx_xy({ix, iy});
                                    size_t indx_ixiy = indx_xy({iy, ix});

                                    return std::real(Q_xy[indx_ixiy]);
                                };

                                double norm_fst = 0.0;
                                double norm_mst = 0.0;

                                integrator_xy->trapz([&density, &Q_xy_f](const size_t &iy, const size_t &ix)
                                                     { return density(iy, ix, Q_xy_f); }, norm_fst);

                                integrator_xy->trapz([&density, &Q_xy_m](const size_t &iy, const size_t &ix)
                                                     { return density(iy, ix, Q_xy_m); }, norm_mst);

                                // std::cout << norm_mst << std::endl;

                                /////////////////////

                                // Fourier

                                /////////////////////

                                std::vector<complex_t> Q_xy_fm(N_xy);

                                for (size_t ix = 0; ix < Nx_qe; ix++)
                                {
                                    for (size_t iy = 0; iy < Ny_qe; iy++)
                                    {
                                        std::complex<double> result = 0.0;
                                        integrator_z->trapz(
                                            [ix, iy, &psi_fst_k, &psi_mst_k, Nx_qe, Ny_qe, norm_fst, norm_mst](const size_t &iz)
                                            {
                                                size_t index = iz * (Ny_qe * Nx_qe) + iy * Nx_qe + ix;
                                                return std::conj(psi_fst_k[index]) * psi_mst_k[index] / (sqrt(norm_fst * norm_mst));
                                            },
                                            result);

                                        // Do something with the result, e.g., store it in a 2D array
                                        size_t indx_ixiy = indx_xy({iy, ix});
                                        Q_xy_fm[indx_ixiy] = result;

                                        // std::cout << "test" << std::endl;
                                    }
                                }

                                std::complex<double> F_S_fm = 0.;

                                auto Fourier_transform = [m, n, &xygrid, &indx_xy, &params](const size_t &iy, const size_t &ix, const std::vector<complex_t> &Q_xy)
                                {
                                    double x = (*xygrid)(ix, iy)[0];
                                    double y = (*xygrid)(ix, iy)[1];

                                    // size_t indx_ixiy = indx_xy({ix, iy});
                                    size_t indx_ixiy = indx_xy({iy, ix});

                                    double Sr = 2. * M_PI / params.a * (1. / sqrt(3.) * (m + n) * x + (m - n) * y);
                                    std::complex<double> PW = exp(I * Sr);

                                    return Q_xy[indx_ixiy] * PW;
                                };

                                integrator_xy->trapz([&Fourier_transform, &Q_xy_fm](const size_t &iy, const size_t &ix)
                                                     { return Fourier_transform(iy, ix, Q_xy_fm); }, F_S_fm);

                                outfile << "band f = " << fst << ", band m = " << mst << ": (" << std::real(F_S_fm) << "," << std::imag(F_S_fm) << ")\n";

                                ////////////////////////////////
                            }
                        }
                        // end of mst and fst loop
                    }
                }

                std::cout << "k point " << ik << std::endl;
            }
        }

        ///////////
    }
    catch (std::string er)
    {
        std::cout << ' ' << er << std::endl;
        std::cout << " Task not accomplished.\n";
        return 1;
    }
    std::cout << "\n Tasks accomplished.\n";
    return 0;
}