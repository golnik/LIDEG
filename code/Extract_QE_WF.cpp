#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>

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

namespace fs = std::filesystem;

// Define file paths
Parameters params;
Parser parser(params);

const std::string QE_path=params.QE_path;
const std::string filepath = QE_path + "wfck2r.oct";
const std::string output_dir = QE_path + "wfc";

std::mutex io_mutex; // Mutex for safe I/O operations
std::mutex queue_mutex; // Mutex for task queue
std::condition_variable cv; // Condition variable for thread synchronization

// Task queue for saving tasks
std::queue<std::tuple<std::vector<std::complex<double>>, int, int>> task_queue;
bool done = false; // Flag to signal thread termination

// Function to parse a line as a complex number
std::complex<double> parseComplex(const std::string &line)
{
    double re, im;
    sscanf(line.c_str(), "(%lf,%lf)", &re, &im);
    return std::complex<double>(re, im);
}

// Function to save a wavefunction to a file
void saveWavefunction(const std::vector<std::complex<double>> &wavefunction, int current_n, int current_k)
{
    std::string filename = output_dir + "/wfc_" + std::to_string(current_n + 1) + "_" + std::to_string(current_k + 1) + ".dat";
    std::ofstream output_file(filename);

    if (output_file.is_open())
    {
        for (const auto &val : wavefunction)
        {
            output_file << val.real() << " " << val.imag() << "\n";
        }
        output_file.close();

        std::lock_guard<std::mutex> lock(io_mutex);
        std::cout << "Saved wavefunction n=" << current_n + 1 << " k=" << current_k + 1 
                  << " to " << filename << " with length " << wavefunction.size() << std::endl;
    }
    else
    {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
    }
}

// Worker function to process tasks from the queue
void worker()
{
    while (true)
    {
        std::tuple<std::vector<std::complex<double>>, int, int> task;
        
        // Retrieve a task from the queue
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            cv.wait(lock, [] { return !task_queue.empty() || done; });
            
            if (done && task_queue.empty())
                return;

            task = std::move(task_queue.front());
            task_queue.pop();
        }
        
        // Process the task
        auto &[wavefunction, current_n, current_k] = task;
        saveWavefunction(wavefunction, current_n, current_k);
    }
}

int main()
{
    // Check if the file exists
    if (!fs::exists(filepath))
    {
        std::cerr << "File not found: " << filepath << std::endl;
        return 1;
    }

    std::ifstream file(filepath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return 1;
    }

    bool isReadingUnkr = false;
    std::vector<int> dims;
    int nr1x, nr2x, nr3x, nbands, nkpoints, wavefunction_size;
    int wavefunction_counter = 0;
    std::string line;

    // First pass: Read dimensions
    while (std::getline(file, line))
    {
        line = line.substr(0, line.find_last_not_of(" \t\n\r\f\v") + 1);

        if (line.find("# name: unkr") != std::string::npos)
        {
            isReadingUnkr = true;
            continue;
        }

        if (isReadingUnkr && line.find("ndims:") != std::string::npos)
        {
            std::getline(file, line);
            std::istringstream dim_stream(line);
            int dim;
            while (dim_stream >> dim)
            {
                dims.push_back(dim);
            }
            if (dims.size() == 5)
            {
                nr1x = dims[0];
                nr2x = dims[1];
                nr3x = dims[2];
                nbands = dims[3];
                nkpoints = dims[4];
                wavefunction_size = nr1x * nr2x * nr3x;
            }
            break;
        }
    }
    
    std::cout << "nr1x: " << nr1x << "\nnr2x: " << nr2x << "\nnr3x: " << nr3x << "\nnbands: " << nbands << "\nnkpoints: " << nkpoints << std::endl;

    fs::create_directories(output_dir);
    std::vector<std::complex<double>> current_wavefunction;

    // Start a pool of worker threads
    int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i)
    {
        threads.emplace_back(worker);
    }

    // Second pass: Read and save each wavefunction incrementally
    while (std::getline(file, line))
    {
        line = line.substr(0, line.find_last_not_of(" \t\n\r\f\v") + 1);

        if (isReadingUnkr && line[0] == '(')
        {
            current_wavefunction.push_back(parseComplex(line));

            if (current_wavefunction.size() == wavefunction_size)
            {
                int current_k = wavefunction_counter / nbands;
                int current_n = wavefunction_counter % nbands;
                wavefunction_counter++;

                // Add the task to the queue
                {
                    std::lock_guard<std::mutex> lock(queue_mutex);
                    task_queue.emplace(current_wavefunction, current_n, current_k);
                }
                cv.notify_one();

                current_wavefunction.clear();
            }
        }
    }

    file.close();

    // Signal worker threads to stop and join them
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        done = true;
    }
    cv.notify_all();

    for (auto &t : threads)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    std::cout << "All wavefunctions have been saved." << std::endl;

    return 0;
}
