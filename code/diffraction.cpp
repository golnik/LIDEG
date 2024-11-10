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
        int tstep = std::stoi(argv[2]) - 1;

        Parameters params;
        Parser parser(params);

        parser.analyze(fname);

        params.print(std::cout);

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
        double time = tgrid[tstep];

        // get vector potential at the given time step
        double Ax = (*Afield_x)(time);
        double Ay = (*Afield_y)(time);

        std::cout << "Ax: " << Ax << " Ay: " << Ay << std::endl;

        /// prepare Black spots
        typedef std::vector<std::vector<double>> BZ_t;
        BZ_t BZ0{{0, 0}};
        BZ_t BZ1{{1, 0}, {1, 1}, {0, 1}, {-1, 0}, {-1, -1}, {0, -1}};
        BZ_t BZ2{{2, 1}, {1, 2}, {-1, 1}, {-2, -1}, {-1, -2}, {1, -1}};
        BZ_t BZ3{{2, 0}, {2, 2}, {0, 2}, {-2, 0}, {-2, -2}, {0, -2}};

        BZ_t BZt{{1, 0}, {-1, 0}};

        std::vector<BZ_t> zones{BZ1};

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

        HexagonalTBModel *tb = new HexagonalTBModel(params.a);

        // create kgrid
        Grid2D *kxygrid;
        Grid2D *Akxygrid;
        double dkx = 0;
        if (params.kgrid_type == kgrid_types::ucell)
        {
            double Ox = 0.;
            double Oy = 0.;
            double b1x = 2. * M_PI / (sqrt(3.) * params.a);
            double b1y = 2. * M_PI / (params.a);
            double b2x = b1x;
            double b2y = -b1y;

            kxygrid = new UCellGrid2D(Ox, Oy, b1x, b1y, params.Nkx, b2x, b2y, params.Nky);
            Akxygrid = new UCellGrid2D(Ox + Ax, Oy + Ay, b1x, b1y, params.Nkx, b2x, b2y, params.Nky);
        }
        else
        {
            throw std::string("Integration is possible only for kgrid_type=ucell!");
        }

        Integrator2D *integrator_kxky = new Integrator2D(kxygrid);
        double SBZ = pow(2. * M_PI, 2.) / (0.5 * sqrt(3.) * params.a * params.a);

        // create graphene model
        Graphene *gm;
        if (params.model == models::hommelhoff)
        {
            double e2p = params.e2p[0];
            double gamma = params.gamma[0];
            double s = params.s[0];
            gm = new GrapheneModel(params.a, e2p, gamma, s, params.Td, E0, E0);
        }
        else if (params.model == models::nlayer)
        {
            gm = new NGraphene(tb, params.nlayers,
                               kxygrid,
                               params.e2p, params.gamma, params.s, params.Td,
                               E0, E0);
        }

        // create rgrid
        Grid2D *xygrid;
        if (params.rgrid_type == rgrid_types::rectan)
        {
            Grid1D *xgrid = new RegularGrid1D(params.xmin, params.xmax, params.Nx);
            Grid1D *ygrid = new RegularGrid1D(params.ymin, params.ymax, params.Ny);

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

            xygrid = new UCellGrid2D(Ox, Oy, a1x, a1y, params.Nx, a2x, a2y, params.Ny);
        }

        Grid1D *zgrid = new RegularGrid1D(params.zmin, params.zmax, params.Nz);

        Integrator2D *integrator_xy = new Integrator2D(xygrid);
        Integrator1D *integrator_z = new Integrator1D(zgrid);

        size_t Nst = gm->nstates();

        /////////////////////////////////////

        MultiIndex indx_xy({params.Nx, params.Ny});
        size_t N_xy = indx_xy.size();

        MultiIndex indx_kxkyst({params.Nkx, params.Nky, Nst});
        size_t N_kxkyst = indx_kxkyst.size();

        /////////////////////////////////////////

        // compute eigenstates in reciprocal space
        std::vector<vector_t> vecs(N_kxkyst);

        for (size_t ist = 0; ist < Nst; ist++)
        {
            for (size_t ikx = 0; ikx < params.Nkx; ikx++)
            {
                for (size_t iky = 0; iky < params.Nky; iky++)
                {
                    // we compute vectors using A-shifted grid!
                    double kx = (*Akxygrid)(ikx, iky)[0];
                    double ky = (*Akxygrid)(ikx, iky)[1];

                    size_t indx_ikxikyist = indx_kxkyst({ikx, iky, ist});
                    vecs[indx_ikxikyist] = gm->get_state(kx, ky, ist);
                }
            }
        }

        // create graphene material
        Orbital *pz = new Pzorb_normal(params.Z);
        Material graphene;

        // generate graphene
        for (size_t il = 0; il < params.nlayers; il++)
        {
            double x = 0.; // shift of the layer depending on stacking

            stacking ABC = params.layers[il];
            switch (ABC)
            {
            case stacking::A:
                x = 0.;
                break;
            case stacking::B:
                x = params.a / sqrt(3.);
                break;
            case stacking::C:
                x = 2. * params.a / sqrt(3.);
                break;
            }

            double z = il * params.d; // position of the layer in z coordinate

            // A..B atoms in graphene layer
            AtomsSet setA = GenerateGraphenePattern(pz, params.a, params.Nclx, params.Ncly, x, 0., z);
            AtomsSet setB = GenerateGraphenePattern(pz, params.a, params.Nclx, params.Ncly, x + params.a / sqrt(3.), 0., z);

            // should be the states with the A field shift
            setA.compute_on_grid(Akxygrid);
            setB.compute_on_grid(Akxygrid);

            graphene.add_atomsset(setA);
            graphene.add_atomsset(setB);
        }

        MultiIndex indx_xyz({params.Nx, params.Ny, params.Nz});
        size_t N_xyz = indx_xyz.size();
        std::vector<complex_t> Q_mf(N_xyz);
        std::vector<complex_t> Q_fn(N_xyz);

        // output diffraction data
        std::string diff_t_fname = params.diffile_fname;
        replace(diff_t_fname, "%it", boost::str(boost::format("%06d") % (tstep + 1)));
        std::ofstream diff_t_out(diff_t_fname);

        std::cout << "Diffraction will be written to: " << diff_t_fname << std::endl;

        diff_t_out << std::scientific;
        diff_t_out << std::setprecision(8);

        // read kspace data
        std::vector<matrix_t> dens_data(Nst);
        std::vector<matrix_t> coh_re_data(Nst * (Nst - 1) / 2);
        std::vector<matrix_t> coh_im_data(Nst * (Nst - 1) / 2);

        std::string dens_t_fname = params.densfile_fname;
        replace(dens_t_fname, "%it", boost::str(boost::format("%06d") % (tstep + 1)));

        size_t col = 0;
        for (size_t ist = 0; ist < Nst; ist++)
        {
            dens_data[ist] = read_2D_from_file<matrix_t>(dens_t_fname, col, params.Nkx, params.Nky);
            col++;
        }

        size_t indx = 0;
        for (size_t ist = 0; ist < Nst; ist++)
        {
            for (size_t jst = ist + 1; jst < Nst; jst++)
            {
                coh_re_data[indx] = read_2D_from_file<matrix_t>(dens_t_fname, col, params.Nkx, params.Nky);
                col++;
                coh_im_data[indx] = read_2D_from_file<matrix_t>(dens_t_fname, col, params.Nkx, params.Nky);
                col++;
                indx++;
            }
        }

        /////////

        diff_t_out << std::setw(25) << time * au2fs;

        //////////

        // The 1st col is time (indx 0)

        // For each Bragg spot start from 2nd col (indx 1)

        // indx [3n + 1] is intra
        // indx [3n + 2] is inter
        // indx [3n + 3] is total
        // where n is nature number

        int t_indx = 0;
        size_t ispot = 0;
        for (size_t izone = 0; izone < nzones; izone++)
        {
            for (auto spot : zones[izone])
            {
                auto m = spot[0];
                auto n = spot[1];

                auto func = [&dens_data, &coh_re_data, &coh_im_data, &params, m, n, nspots, nzones, &zones](const size_t &ikx, const size_t &iky)
                {
                    std::vector<std::complex<double>> res_k(3, 0.);
                    std::complex<double> F_S[8][8];

                    // Initialize all elements to zero
                    for (int i = 0; i < 8; ++i)
                    {
                        for (int j = 0; j < 8; ++j)
                        {
                            F_S[i][j] = std::complex<double>(0.0, 0.0);
                        }
                    }

                    // Construct the file path and open the file
                    std::string file_path = "/xdisk/ngolubev/mingruiyuan/QE_diffr/transition_re/" + std::to_string(ikx * params.Nkx + iky + 1) + "_rearrange.dat";
                    std::ifstream infile(file_path);

                    if (!infile.is_open())
                    {
                        std::cerr << "File " << file_path << " not found!" << std::endl;
                    }

                    // Create the regex pattern with m and n values
                    std::ostringstream pattern_stream;
                    pattern_stream << "S " << m << "," << n << R"( band f = (\d+), band m = (\d+): \((-?[\d\.eE+-]+),(-?[\d\.eE+-]+)\))";
                    std::string pattern = pattern_stream.str();
                    std::regex sum_pattern(pattern);
                    std::smatch matches;
                    std::string line;

                    // Read the file line by line
                    while (std::getline(infile, line))
                    {
                        if (std::regex_search(line, matches, sum_pattern))
                        {
                            int band_f = std::stoi(matches[1]) - 1;
                            int band_m = std::stoi(matches[2]) - 1;

                            // Ensure band indices are within bounds
                            if (band_f >= 0 && band_f < 8 && band_m >= 0 && band_m < 8)
                            {
                                double real_part = std::stod(matches[3]);
                                double imag_part = std::stod(matches[4]);

                                // Store the complex number in F_S
                                F_S[band_f][band_m] = std::complex<double>(real_part, imag_part);
                            }
                            else
                            {
                                std::cerr << "Band index out of bounds: band f = " << band_f + 1 << ", band m = " << band_m + 1 << std::endl;
                            }
                        }
                    }

                    infile.close();

                    // Continue with calculations as before
                    for (size_t mst = 3; mst < 5; mst++)
                    {
                        for (size_t fst = 3; fst < 5; fst++)
                        {
                            for (size_t nst = 3; nst < 5; nst++)
                            {
                                std::complex<double> rho[2][2];
                                rho[0][0] = dens_data[0](ikx, iky);
                                rho[1][1] = dens_data[1](ikx, iky);
                                rho[0][1] = coh_re_data[0](ikx, iky) + I * coh_im_data[0](ikx, iky);
                                rho[1][0] = coh_re_data[0](ikx, iky) - I * coh_im_data[0](ikx, iky);

                                if (mst == nst) // intra band
                                {
                                    res_k[0] += rho[nst - 3][mst - 3] * std::conj(F_S[fst][mst]) * F_S[fst][nst];
                                }

                                if (mst != nst) // inter band
                                {
                                    res_k[1] += rho[nst - 3][mst - 3] * std::conj(F_S[fst][mst]) * F_S[fst][nst];
                                }

                                res_k[2] += rho[nst - 3][mst - 3] * std::conj(F_S[fst][mst]) * F_S[fst][nst];
                            }
                        }
                    }

                    vector<double> res(3, 0.);
                    for (size_t i = 0; i < res_k.size(); i++)
                    {
                        res[i] = std::real(res_k[i]);
                    }

                    return res;
                };

                vector<double> res(3, 0.);

                integrator_kxky->trapz(func, res);

                res *= 2. / SBZ;

                ispot++;

                std::cout << "Ispot         : " << ispot << std::endl;

                diff_t_out << std::setw(25) << res(0);

                diff_t_out << std::setw(25) << res(1);

                diff_t_out << std::setw(25) << res(2);
            }
        }

        ///////////

        diff_t_out.close();
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