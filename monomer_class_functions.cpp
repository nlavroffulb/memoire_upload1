#include "monomer_class.h"

monomer::monomer(std::vector<double> position, int identity, std::string base_type) {
    P = position;
    id = identity;
    base = base_type;
}

std::vector<double> monomer::get_position() {
    return P;
}
bool monomer::monomer_is_linked() {
    if (linked_monomer == nullptr) {
        std::cout << "False" << "\n";
        linked = false;
        return false;
    }
    else {
        std::cout << "True" << "\n";
        linked = true;
        return true;
    }
}


void monomer::link_to(monomer* link_to_or_null) {
    //what we're doing to this monomer we have to do to the monomer being linked

    if (link_to_or_null == nullptr) {
        linked = false;
        linked_monomer = nullptr;
    }
    else {
        linked = true;
        linked_monomer = link_to_or_null;

        link_to_or_null->linked = true;
        link_to_or_null->linked_monomer = this;

    }


}


void monomer::change_position(std::vector<double> new_position)
{
    P = new_position;
}