#include <iostream>
#include<vector>
#include "polymer_class.h"
#include "polymer_generation.h"
#include<time.h>//set seed.
#include <fstream>
#include <chrono>//measure performance time


int main()
{
    srand(time(NULL));

    //endpoint definitions
    std::vector<double> initial_position{0,0,0.1};
    std::vector<double> N_position{ 0,0,30 };

    //simple example of growing the polymer once and can output txt file to view it in ovito. can also measure
    // performance.
    // 
    //auto start = std::chrono::high_resolution_clock::now();
    //std::vector<std::vector<double>> p{ grow_chain(N_position,initial_position,70) };
    //output_for_ovito(p);
    //auto stop = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    //std::cout << duration.count() << std::endl;

    //********************************************************************************
    //********************************************************************************
    //********************************************************************************

    //code to investigate e2e distance of middle monomer for polymer of length 71. 
    // results outputted to text file.
    // 
    //auto start = std::chrono::high_resolution_clock::now();
    //static const int N_sims{ 10000 };
    //std::vector<double> sim_results(N_sims);
    //auto stop = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    //std::cout << duration.count() << std::endl;
    //for (int i{ 0 }; i < N_sims; i++) {
    //    std::cout << i << std::endl;
    //    std::vector<std::vector<double>> p{ grow_chain(N_position,initial_position,70) };
    //    //sim_results.push_back(dist_2_points3d(p[59], initial_position));
    //    //sim_results.push_back(dist_2_points3d(p[59], N_position));
    //    sim_results[i]=dist_2_points3d(p[35], initial_position);
    //    //std::cout << dist_2_points3d(p[25], initial_position) << std::endl;
    //    //sim_results[elements-1-i] = dist_2_points3d(p[25], N_position);
    //    //std::cout << dist_2_points3d(p[25], N_position) << std::endl;
    //}
    //
    //std::ofstream vector("sim_results.txt");
    //for (int k{ 0 }; k < N_sims; k++) {
    //    vector << sim_results[k] << std::endl;;
    //********************************************************************************
    //********************************************************************************
    //********************************************************************************

    //using the polymer object and testing some of the functions.
    // 
    //polymer RNA(p);
    //std::cout << RNA.regions_compatible("UCAGA", "UCUGA") << std::endl;
    //std::vector < std::vector<std::vector<int>>> s{ RNA.structure_search() };
    //RNA.sample_structure(s);
    //RNA.structure_search();
    //output_for_ovito(p);

}
