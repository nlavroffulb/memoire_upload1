#pragma once
#pragma once

#include "monomer_class.h"
//#include "structure_class.h"
//#include "math_functions.h"



class polymer {
protected:
    monomers chain;

private:
    int N{ 0 };

public:

    //constructors
    polymer() = default;
    polymer(int num);

    polymer(std::vector<double> beginning_position, std::vector<double> end_position);
    polymer(std::vector<std::vector<double>> generated_positions);

    //general functions of the polymer and given monomer in the chain
    int chain_length();
    void get_monomer_position(int id);
    void add_monomer(std::vector<double> position, std::string base_type);
    std::vector<std::vector<double>> generate_vector_of_monomer_positions();
    void print_monomer_positions();
    void output_for_ovito();

    //operator functions
    monomer* operator[](int i);

    std::string random_base();


    //compatibility
    bool base_pairing_rules(char b1, char b2);
    std::vector<int> position_of_region(int i, int r);
    bool regions_compatible(std::string r1, std::string r2);

    std::vector<std::vector<std::vector<int>>> structure_search();
    std::vector<int> sample_structure(std::vector<std::vector<std::vector<int>>> potential_structures);

    bool compatible_bases(monomer* monomer_i, monomer* monomer_j);
    void create_structure_from_compatible_region(std::vector<int> compatible_region);
    ////////////////////////////////////////////////////////////////////////////////////////

    //re-generating polymer with structure
    std::vector<double> rotate_helix(std::vector<double> u, std::vector<double> v, double the);
    void generate(std::vector<int> compatible_monomers, std::vector<double> v, std::vector<double> u, std::vector<double> x);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void classify_compatible_regions(std::vector<std::vector<int>> search_results);
    void link_move();

    //destructor
    ~polymer() { std::cout << "Polymer destructor called.\n"; }



};

