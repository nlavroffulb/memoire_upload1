#include "polymer_generation.h"

// l = number of segments to be regrown
// i ranges from 1 to l
// the total number of monomers is l+1 (counting the start and end monomers)

double const a{ 1.0 };


double lennard_jones_ij(std::vector<double> position_i, std::vector<double> position_j) {
	if (!position_i.empty() && !position_j.empty()) {
		double r_ij{ dist_2_points3d(position_i,position_j) };

		if (r_ij <= a * pow(2, 1.0 / 6.0)) {
			return 4 * (pow(a / r_ij, 12) - pow(a / r_ij, 6)) + 1;
		}
		else if (r_ij > a * pow(2, 1.0 / 6.0)) {
			return 0;
		}

	}
}

//excluded volume potential at a given position (ie the position of a monomer)
double u_r(std::vector<double> r, std::vector<std::vector<double>> positions_ij) {
	double summand{ 0 };
	for (int j{ 0 }; j < positions_ij.size(); j++) {
		//sometimes we will calculate u(r) and positions_ij will be basically empty.
		//issue is that we declare its size at the start, so it can still loop through the elements
		//even if the elements themselves are empty.

		if (positions_ij[j] == r) {
			continue;
		}
		else {
			summand += lennard_jones_ij(r, positions_ij[j]);
		}
	}
	return summand;
}

//rosenbluth sampling selection probability for a given trial position. the argument assumes that the energy
// of the trial position was calculated before this function is called.
double selection_prob(std::vector<double>& trial_energies, double& trial_energy) {
	double normalization{0};
	for (auto i : trial_energies) {
		normalization += exp(-i);
	}
	return exp(-trial_energy) / normalization;
}

// function which takes as argument all the trial positions that have been generated
//and selects one of them using a roulette method.
std::vector<double> select_trial(std::vector<std::vector<double>>&trial_positions, std::vector<std::vector<double>> &existing_positions)
{

	std::vector<double> bins(trial_positions.size() + 1);

	bins[0] = 0;
	double summand{ 0 };
	for (int j{ 0 }; j < trial_positions.size(); j++) {
		summand += exp(-u_r(trial_positions[j],existing_positions));
		bins[j + 1] = summand;
	}
	double rn{ rand2(0,1)*summand };
	int selection{ 0 };

	for (int k{ 1 }; k < bins.size(); k++) {
		if (rn < bins[k]) {
			selection = k - 1;
			break;//once we've found the corresponding bin we can break the loop and return the selected trial.
		}
	}
	return trial_positions[selection];
}


double ideal_chain_pdf(double e2e_distance, int n_segments)//'yamakawa' function.
// protected against factorial overflow if n>175.
{
	if (n_segments < 170) {
		double sum_value{ 0 };
		for (int k{ 0 }; k <= (n_segments - (e2e_distance / a)) / 2; k++) {
			sum_value += pow(-1, k) * choose(n_segments, k) * pow(n_segments - 2.0 * k - (e2e_distance / a), n_segments - 2);
		}
		return sum_value / (pow(2, n_segments + 1) * factorial(n_segments - 2) * pi * pow(a, 2) * e2e_distance);

	}
	else if (n_segments > 170) {//needs to be updated. will have a lookup table to a gaussian distribution.
		return pow(1.5 * pi * n_segments, 3.0 / 2.0) * exp(-3.0 * pow(e2e_distance, 2) / (2 * n_segments));
	}
}

//this is the main probability distribution that is called when the polymer is grown. uses ideal_chain_pdf
double segment_grow_prob(double new_e2e, double old_e2e, int l_minus_i) {
	//return ideal_chain_pdf(new_e2e, l_minus_i) / (4*pi*ideal_chain_pdf(old_e2e, l_minus_i +1));
	return ideal_chain_pdf(new_e2e, l_minus_i);

}

//check for overstretch since all the segments have a fixed length.
bool overstretch(int l_minus_i, double trial_e2e) {
	if (trial_e2e > (l_minus_i)* a) {
		return true;
	}
	else { return false; }
}

//needs to be updated with ran2 recipe. 
//random number generator
double rand2(double min_value, double max_value) {

	return double(((max_value - min_value) * rand()) / RAND_MAX) + min_value;
}

//generate a trial jump from the current position where we trial a position for the next monomer to be added
void sample_jump_direction(std::vector<double> &jump, double bond_distance) {
	bond_distance = a;

	double phi{ rand2(0,2 * pi) }, theta{ acos(1-2*rand2(0,1)) };//correct sampling of cos theta using 
	//inverse transform sampling

	jump = { a * sin(theta) * cos(phi),a * sin(theta) * sin(phi),a * cos(theta) };

}


std::vector<double> rejection_sample(std::vector<double> initial, std::vector<double> N_position, int segments_to_regrow) {

	//this if statement might not be necessary because we try to be careful in what follows
	//to treat crankshaft insertion differently and never call this function once we reach the last
	//monomer to be grown.
	if (segments_to_regrow > 1) {
		//pdf max. needed for von Neumann rejection sampling
		double init_e2e{ vector_modulus(vector_subtraction(N_position,initial)) };//need for pdf:max
		double pdf_max{ segment_grow_prob(init_e2e - a,init_e2e,segments_to_regrow) };
		double Y{ rand2(0,pdf_max) };//for von Neumann rejection sampling

		//trial jump
		std::vector<double> jump(3);
		sample_jump_direction(jump, a);
		std::vector<double> trial_position{ vector_addition(initial,jump) };
		double trial_e2e{ vector_modulus(vector_subtraction(N_position,trial_position)) };


		int over_count{ 0 };//variable to keep track of the number of rejections

		//we reject based on overstretch condition and based on von Neumann.
		while (overstretch(segments_to_regrow, trial_e2e) == true || Y > segment_grow_prob(trial_e2e, init_e2e, segments_to_regrow)) {
			//std::cout << "rejected overstretch" << std::endl;
			Y = rand2(0, pdf_max);

			//generate new trial positions while trials are rejected
			sample_jump_direction(jump, a);
			trial_position = vector_addition(initial, jump);
			trial_e2e = vector_modulus(vector_subtraction(N_position, trial_position));
			over_count++;
			if (over_count > 10000) {
				std::cout << "stuck" << std::endl;
			}
		}

		return trial_position;

	}
	else {
		std::cout << "The last segment needs to be regrown by crankshaft insertion" << std::endl;
		
	}

}

// method for crankshaft insertion... view latex.
//there are a lot of variable definitions... will be made more compact at some point.
std::vector<double> crankshaft_insertion(std::vector<double> start_pos, std::vector<double> N_pos) {

	std::vector<double> r_e2e{ vector_subtraction(N_pos,start_pos) };//the e2e distance
	//between the N-2 monomer and the Nth monomer (the endpoint).


	double mod_r_e2e{ vector_modulus(r_e2e) };//distance between N-2 and N

	double r_c{ sqrt(a * a - pow(mod_r_e2e / 2,2)) };//r_c is the radius of the circle at the top of the cone

	//we sample x and y coordinates on the circle of the cone. specifically we sample theta of the circle
	// in polar coordinates.
	// 
	// then use definition of dot product to solve for z. 
	//Then vector addition/subtraction to get back to our coordinate system

	//sample x and y
	r_e2e = normalize(r_e2e);
	double theta = rand2(0, 2 * pi);
	double phi = acos(mod_r_e2e / (2 * a));
	double x{ cos(theta) }, y{ sin(theta) };

	//solve for z
	double z{ -1 / r_e2e[2] * (x * r_e2e[0] + y * r_e2e[1]) };

	std::vector<double> v_c{ x,y,z };
	v_c = normalize(v_c);
	r_e2e = multiplication_by_scalar(mod_r_e2e / 2, r_e2e);
	v_c = multiplication_by_scalar(r_c, v_c);
	//vector addition to get back to our coordinate system
	std::vector<double> r_Nminus1{ vector_addition(r_e2e,v_c) };

	r_Nminus1 = vector_addition(start_pos, r_Nminus1);

	return r_Nminus1;
};

// main function to generate a polymer with fixed endpoints and constant segment length.
//note the final vector of positions includes the start and endpoints. 
//so for n segments to grow there are n+1 elements in the final vector of positions that we generate.
std::vector<std::vector<double>> grow_chain(std::vector<double> starting_end, std::vector<double> ending_end, int l) {
	std::vector<std::vector<double>> positions(l+1);
	std::vector<double> generated_position(3);
	std::vector<double> old_position{ starting_end };
	positions[0]=starting_end;
	int i_segments;

	//i is the monomer that is being grown. 
	// l is the total number of segments to be regrown.
	// if you have 4 segments in total for example then you need
	//to add 3 monomers in addition to the endpoints so i=1,2,3 (l-1 = 3). 
	// Note that including the endpoints we have 5 monomers in total in the flexible chain in the example. 
	for (int i=1; i <= l-1; i++) {
		//std::cout <<"regrowth stage " << i << std::endl;
		if (l - i > 1) {
			int trials{ 100 };
			std::vector<std::vector<double>> trial_positions(trials);
			std::vector<double> trial_position(3);
			i_segments = l - i;

			//rosenbluth sampling excluded volume
			for (int k{ 0 }; k < trials; k++) {
				//std::cout << k << std::endl;
				trial_position = rejection_sample(old_position, ending_end, i_segments);
				trial_positions[k] = trial_position;
			}
			generated_position = select_trial(trial_positions,positions);
			old_position = generated_position;

			positions[i] = generated_position;

		}
		else if(l-i==1) {

			//rosenbluth sampling not yet implemented for crankshaft insertion.
			generated_position = crankshaft_insertion(positions[l-2], ending_end);
			positions[l-1] = generated_position;

		}
	}
	positions[l] = ending_end;

	//uncomment to print out the generated positions
	//for (int i = 0; i < positions.size(); i++)
	//{
	//	std::cout << i;
	//	std::cout << "{";
	//	for (int j = 0; j < positions[i].size(); j++)
	//	{
	//		std::cout << positions[i][j] << " ";
	//	}
	//	std::cout << "}\n";
	//}

	return positions;
}

// another function to output for ovito. Defined again in the polymer class. Could delete one or the other.
// some of the interesting regions are coloured differently.
void output_for_ovito(std::vector<std::vector<double>> positions) {

	std::string file_name{ "polymer " + std::to_string(rand2(0,1)) + ".txt" };
	std::ofstream ovito_output(file_name);
	
	ovito_output << positions.size() << "\n\n";
	for (int i{ 0 }; i < positions.size(); i++) {
		std::string line{ "" };
		if (i == 0) {//startpoint monomer
			line += "He";
		}
		else if (i == positions.size() - 2) {//crankshaft insertion monomer
			line += "Mg";
		}
		else if (i == positions.size() - 1) {//endpoint
			line += "He";
		}

		else if (i % 2 == 0) {
			line += "C";
		}
		else if (i % 2 != 0) {
			line += "O";
		}
		for (int j{ 0 }; j < positions[i].size(); j++) {
			line += " " + std::to_string(positions[i][j]);
		}
		if (i == positions.size() - 1) {
			ovito_output << line << std::endl;;
		}
		else {
			ovito_output << line << "\n";
		}
	}

	ovito_output.close();
}


