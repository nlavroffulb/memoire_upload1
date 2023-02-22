#pragma once

#include <vector>
#include <iostream>

class monomer {
private:
    std::vector<double> P{ 0,0,0 };
    bool linked{ false };
    monomer* linked_monomer{ nullptr };
    std::string base{ "" };

    //structure* part_of_structure{ nullptr };
public:
    int id{ 0 };

    //constructors
    monomer() = default;
    monomer(std::vector<double> position, int identity, std::string base_type);

    // access member variables
    std::vector<double> get_position();
    bool monomer_is_linked();

    int get_id() { return id; }
    std::string get_base() { return base; }

    void link_to(monomer* link_to_or_null);

    void change_position(std::vector<double> new_position);
    //destructor
    ~monomer() { std::cout << "Monomer destructor called\n"; }
};

typedef std::vector<monomer*> monomers;
