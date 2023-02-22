#pragma once
#include<vector>
#include "math_functions.h"
#include <fstream>
#include<string>
#include<iostream>


double selection_prob(std::vector<double>& trial_probabilities, double& trial_i);
std::vector<double> select_trial(std::vector<std::vector<double>>& trial_positions, std::vector<std::vector<double>>& existing_positions);

static const double pi{ atan(1) * 4 };
double rand2(double min_value, double max_value);

double ideal_chain_pdf(double e2e_distance, int n_segments);
double segment_grow_prob(double new_e2e, double old_e2e, int l_minus_i);
bool overstretch(int l_minus_i, double new_e2e);

void sample_jump_direction(std::vector<double> &jump, double bond_distance);
std::vector<double> crankshaft_insertion(std::vector<double> initial, std::vector<double> N_position);
std::vector<double> rejection_sample(std::vector<double> initial, std::vector<double> N_position, int number_of_segments);

std::vector<std::vector<double>> grow_chain(std::vector<double> starting_end, std::vector<double> ending_end, int l);
void output_for_ovito(std::vector<std::vector<double>> positions);






