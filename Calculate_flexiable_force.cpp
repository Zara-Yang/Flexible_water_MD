#include <iostream>
#include <string>
#include <time.h> 
#include <math.h>   
#include <stdio.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

#define pi 3.141592653589793
#define sqr_pi 1.772453850905516
#define k_l	229074					// unit : 10 J/mol/A**2
#define k_theta 20878.16			// unit : 10 J.mol/rad**2
#define standard_l 1				// unit : A
#define standard_theta 1.91061217	// unit : rad

#define ep0p12 8.8541878128
#define kbp23 1.380649
#define elechp19 1.602176634
#define nan23 6.02214076
#define alpha_eel 0.2222222
#define r4pie_eel 138935.45764438206
#define sig_lj 3.166
#define eps_lj 65
#define b1_sh 6
#define b2_sh 10
#define UnitScale 10
#define nAtom 3072

#define ChargeO -0.8476
#define ChargeH +0.4238
#define MassO 15.9994
#define MassH 1.00794

#define deltaT 0.001
#define TargetTemp 300
#define Totalstep 10000

#define Box_x 31.32019
#define Box_y 31.32019
#define Box_z 31.32019

#define mchk 10

typedef vector<double> VEC3;
typedef vector < vector < double >> MAT3;

VEC3 operator+(const VEC3& v1, const VEC3& v2) {
	VEC3 result(3);
	for (int index = 0; index < 3; index++) {
		result[index] = v1[index] + v2[index];
	}
	return(result);
}

VEC3 operator*(const double c, const VEC3& v) {
	VEC3 result(3);
	for (int index = 0; index < 3; index++) {
		result[index] = c * v[index];
	}
	return(result);
}

VEC3 operator*(const VEC3& v1, const VEC3& v2) {
	VEC3 result(3);
	for (int index = 0; index < 3; index++) {
		result[index] = v1[index] * v2[index];
	}
	return(result);
}

VEC3 operator-(const VEC3& v1, const VEC3& v2) {
	VEC3 result(3);
	for (int index = 0; index < 3; index++) {
		result[index] = v1[index] - v2[index];
	}
	return(result);
}

MAT3 operator*(const double c, const MAT3& m) {
	MAT3 result(3, VEC3(3));
	for (int x = 0; x < 3; x++) {
		for (int y = 0; y < 3; y++) {
			result[x][y] = m[x][y] * c;
		}
	}
	return(result);
}

VEC3 operator*(const MAT3& m, const VEC3& v) {
	VEC3 result(3);
	result[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
	result[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
	result[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
	return(result);
}

MAT3 operator-(const MAT3& n, const MAT3& m) {
	MAT3 result(n.size(), VEC3(n[0].size()));
	for (int x = 0; x < n.size(); x++) {
		for (int y = 0; y < n[0].size(); y++) {
			result[x][y] = n[x][y] - m[x][y];
		}
	}
	return(result);
} 

MAT3 operator+(const MAT3& v1, const MAT3& v2) {
	int size_x = v1.size();
	int size_y = v1[0].size();
	MAT3 result(size_x, VEC3(size_y));

	for (int index_x = 0; index_x < size_x;index_x++){
		for (int index_y = 0; index_y < size_y; index_y++) {
			result[index_x][index_y] = v1[index_x][index_y] + v2[index_x][index_y];
		}
	}
	return(result);
}

ostream& operator<<(ostream& out, const MAT3& v) {
	for (int x = 0; x < v.size(); x++) {
		for (int y = 0; y < 3; y++) {
			out << v[x][y] << " ";
		}
		cout << endl;
	}
	return(out);
}

double dot(VEC3 v1, VEC3 v2) {
	double result = 0;
	for (int index = 0; index < 3; index++) {
		result += v1[index] * v2[index];
	}
	return(result);
}

double sum(VEC3 v){
	double result = 0;
	for (int index = 0; index < v.size(); index ++){
		result += v[index];
	}
	return(result);
}

double norm(VEC3 v){
	double result = 0;
	for (int index = 0; index < v.size(); index ++){
		result += v[index] * v[index];
	}
	result = sqrt(result);
	return(result);
}

void shift_parameter_invr(double alpha_i, double b1_i, double b2_i, double& c0_o, double& c3_o, double& c4_o){
	///=============================================
	/// get shift function parameters c0 c3 and c4 for u(r) = 1 / r^\alpha
	/// ------------------------------------------------------------------
	/// math : 
	///				  |       1
	///				  | ------------- + c0										for r < b1
	/// 			  |    r^\alpha 	 
	///				  |
	///				  |      1
	/// u_shift(r) =  | ------------- + c3 * (r - b1)^3 + c4*(r - b1)^4 + c0	for b1 < r < b2
	///				  |   r^\alpha
	///				  |
	///				  |			
	///               |		0													for b2 < r
	///				  |  
	///	
	/// u_shift potential and force ( negative derivative of u_shift ) is continuous and differentiable. 
	/// =============================================
	double inv_buffer = (b2_i - b1_i) * (b2_i - b1_i) * pow(b2_i, (alpha_i + 2));
	inv_buffer = 1 / inv_buffer;
	c3_o = (1. / 3.) * ((alpha_i + 4.) * b2_i - (alpha_i + 1.) * b1_i) * inv_buffer;
	c3_o = c3_o * alpha_i;
	c4_o = -(1. / 4.) * ((alpha_i + 3.) * b2_i - (alpha_i + 1.) * b1_i) * inv_buffer / (b2_i - b1_i);
	c4_o = c4_o * alpha_i;

	c0_o = -1. / pow(b2_i, alpha_i) - c3_o * pow((b2_i - b1_i), 3) - c4_o * pow((b2_i - b1_i), 4);
}

void shift_parameter_LJ(double sigma, double b1_i, double b2_i, double& c0_o, double& c3_o, double& c4_o) {
	/// ====================================================
	/// shift function of Lennard-Jones potential
	/// 
	/// A = sigma / r
	///			
	///								| A^12 - A^6 + c0										for r < b1
	///								| 
	/// u_shift / (4 * epsilon) =	| A^12 - A^6 + c3 * (r - b1)^3 + c4 * (r - b2)^4 + c0	for b1 < r < b2
	///								|
	///								| 0														for b2 < r
	/// 
	/// ====================================================
	double c0_6, c3_6, c4_6;
	double c0_12, c3_12, c4_12;
	double sigma_6 = pow(sigma, 6);
	double buffer_6 = sigma_6 / pow(b1_i, 6);

	if (abs(b2_i - b1_i) < 1e-10) {
		c0_o = -buffer_6 * buffer_6 + buffer_6;
		c3_o = 0;
		c4_o = 0;
	}
	else {
		shift_parameter_invr(6, b1_i, b2_i, c0_6, c3_6, c4_6);
		shift_parameter_invr(12, b1_i, b2_i, c0_12, c3_12, c4_12);
		c0_o = sigma_6 * (sigma_6 * c0_12 - c0_6);
		c3_o = sigma_6 * (sigma_6 * c3_12 - c3_6);
		c4_o = sigma_6 * (sigma_6 * c4_12 - c4_6);
	}
}

void shift_parameter_erfc(double alpha, double b1_i, double b2_i, double& c0_o, double& c3_o, double& c4_o) {
	/// ==========================================================
	/// Get shift function parameter c0 c3 c4 of u(r) = erfc(alpha * r) / r 
	/// ==========================================================
	/// 
	double sqrpi = 1.772453850905516;
	double erfc_b2 = erfc(alpha * b2_i);
	double exp_b2 = exp(-(alpha * b2_i) * (alpha * b2_i));
	
	double erfc_b2_prime = -erfc_b2 / (b2_i * b2_i) - (2 * alpha) / (b2_i * sqrpi) * exp_b2;
	double erfc_b2_prime2 = 2 * erfc_b2 / pow(b2_i, 3) + 4 * alpha / (b2_i * b2_i * sqrpi) * exp_b2 + 4 * pow(alpha, 3) / sqrpi * exp_b2;

	double b12 = b2_i - b1_i;

	c3_o = (erfc_b2_prime2 / 3 * b12 - erfc_b2_prime) / (b12 * b12);
	c4_o = -(erfc_b2_prime2 * 0.5 * b12 - erfc_b2_prime) / (2 * b12 * b12 * b12);
	c0_o = -erfc_b2 / b2_i - c3_o * pow(b12, 3) - c4_o * pow(b12, 4);
}

VEC3 PeriodicBoundaryCondition(VEC3& coord, VEC3& Box){
	VEC3 result = {0, 0, 0};
	result[0] = coord[0] - round(coord[0]/Box[0])*Box[0];
	result[1] = coord[1] - round(coord[1]/Box[1])*Box[1];
	result[2] = coord[2] - round(coord[2]/Box[2])*Box[2];
	return(result);
}

void Intermolecular_force(double natom, VEC3& charge, MAT3& coord, MAT3& force, VEC3& box_range, double& potential, double shift_b1,double shift_b2,double lj_sigma,double lj_eps,double erfc_alpha, double erfc_r4pie, double lj_c0, double lj_c3, double lj_c4,double erfc_c0,double erfc_c3,double erfc_c4){
    /// ================================================
    /// Functions to compute force and potential for shifted lennard-jones and erfc electrostatics
    ///                      f         p             sh_     l       j_        e    el_
    /// under periodic boundary condition in 3 dimensions H2O system
    ///       p                              3
    /// fp_shljeel_p3 works for H20 spce model only
    ///
    /// Already check with benchmark.
	/// 										2022/1/9
    /// =================================================

	double lj_sigma_sq, lj_eps_4, shift_b1_sq, shift_b2_sq;
	double lj_eps_48, erfc_alpha_sq;
	double r4pie_ec0, r4pie_ec3, r4pie_ec4;
	double eps_lc0, eps_lc3, eps_lc4;

	double qij,r4qij;
	double erfc_r,exp_r;

	lj_sigma_sq = lj_sigma * lj_sigma;
	lj_eps_4 = 4. * lj_eps;
	shift_b1_sq = shift_b1 * shift_b1;
	shift_b2_sq = shift_b2 * shift_b2;

	lj_eps_48 = 12 * lj_eps_4;
	erfc_alpha_sq = erfc_alpha * erfc_alpha;
	r4pie_ec0 = erfc_r4pie * erfc_c0;
	r4pie_ec3 = erfc_r4pie * erfc_c3;
	r4pie_ec4 = erfc_r4pie * erfc_c4;
	eps_lc0 = lj_c0 * lj_eps_4;
	eps_lc3 = lj_c3 * lj_eps_4;
	eps_lc4 = lj_c4 * lj_eps_4;

	for (int index = 0; index < natom; index++){
		force[index] = {0, 0, 0};
	}
	potential = 0;

	VEC3 coord_i(3),coord_j(3),vec_rij(3);
	int atom_type_i, atom_type_j;
	double norm_rij,normsq_rij;
	double uij, wij;
	double sig_unit2,sig_unit6;
	double rb1_diff,rb1_diff_sq;
	double fij;
	double mol1,mol2;
	VEC3 fij_vec;

	for (int index_i = 0; index_i < natom - 1; index_i ++){
		for (int index_j = index_i + 1; index_j < natom; index_j++){
			mol1 = floor(index_i/3.);
			mol2 = floor(index_j/3.);
			if(mol1 == mol2){continue;}
			coord_i = coord[index_i];
			coord_j = coord[index_j];
			vec_rij = coord_i - coord_j;
			
			vec_rij[0] = vec_rij[0] - round(vec_rij[0]/box_range[0])*box_range[0];
			vec_rij[1] = vec_rij[1] - round(vec_rij[1]/box_range[1])*box_range[1];
			vec_rij[2] = vec_rij[2] - round(vec_rij[2]/box_range[2])*box_range[2];

			normsq_rij = dot(vec_rij,vec_rij);
			norm_rij = sqrt(normsq_rij);
			if (normsq_rij < shift_b2_sq){
				qij = charge[index_i] * charge[index_j];
				r4qij = erfc_r4pie * qij;
				erfc_r = erfc(erfc_alpha * norm_rij);
				exp_r = exp(-erfc_alpha_sq * normsq_rij);
				uij = (erfc_r / norm_rij) * r4qij + r4pie_ec0 * qij;
				wij = (erfc_r / norm_rij + 2 * erfc_alpha / sqr_pi * exp_r) * r4qij;
				if ((index_i % 3) == 0 && (index_j % 3) == 0){
					sig_unit2 = lj_sigma_sq / normsq_rij;
					sig_unit6 = sig_unit2 * sig_unit2 * sig_unit2;

					uij += sig_unit6 * (sig_unit6 - 1) * lj_eps_4 + eps_lc0;
					wij += sig_unit6 * (sig_unit6 - 0.5) * lj_eps_48;
				}
				
				if (normsq_rij > shift_b1_sq){
					rb1_diff = norm_rij - shift_b1;
					rb1_diff_sq = rb1_diff * rb1_diff;

					uij = uij + (r4pie_ec3 * rb1_diff + r4pie_ec4 * rb1_diff_sq) * rb1_diff_sq * qij;
					wij = wij - (3 * r4pie_ec3 + 4 * r4pie_ec4 * rb1_diff) * rb1_diff_sq * qij * norm_rij;
					if ((index_i % 3) == 0 && (index_j % 3) == 0){
						uij = uij + (eps_lc3 * rb1_diff + eps_lc4 *rb1_diff_sq) * rb1_diff_sq;
						wij = wij - (3 * eps_lc3 + 4 * eps_lc4 * rb1_diff) * rb1_diff_sq * norm_rij;
					}
				}
				fij = wij / normsq_rij;
				fij_vec = fij * vec_rij;
				force[index_i] = force[index_i] + fij_vec;
				force[index_j] = force[index_j] - fij_vec;
				potential += uij;
			}
		}
	}
}

void Intramolecular_force_for_single(VEC3& Oxyz, VEC3& H1xyz, VEC3& H2xyz, VEC3& Box, VEC3& Force_O, VEC3& Force_H1, VEC3& Force_H2, double& potential){
	potential = 0;
	Force_O  = {0, 0, 0};
	Force_H1 = {0, 0, 0};
	Force_H2 = {0, 0, 0};

	VEC3 vec_OH1 = Oxyz  - H1xyz;
	VEC3 vec_OH2 = Oxyz  - H2xyz;
	VEC3 vec_H12 = H1xyz - H2xyz;

	vec_OH1 = PeriodicBoundaryCondition(vec_OH1, Box);
	vec_OH2 = PeriodicBoundaryCondition(vec_OH2, Box);
	vec_H12 = PeriodicBoundaryCondition(vec_H12, Box);

	double norm_OH1 = norm(vec_OH1);
	double norm_OH2 = norm(vec_OH2);
	double norm_H12 = norm(vec_H12);

	double cos_theta = ( norm_OH1 * norm_OH1 + norm_OH2 * norm_OH2 - norm_H12 * norm_H12 ) / ( 2 * norm_OH1 * norm_OH2 );
	double theta = acos(cos_theta);
	potential = potential + k_l * ( norm_OH1 - standard_l ) * ( norm_OH1 - standard_l );
	potential = potential + k_l * ( norm_OH2 - standard_l ) * ( norm_OH2 - standard_l );
	potential = potential + k_theta * (theta - standard_theta) * (theta - standard_theta);
	double partial_theta_cos = ( -1 ) / sqrt( 1 - cos_theta * cos_theta);

	double partial_cos_OH1 = (   norm_OH1 * norm_OH1 - norm_OH2 * norm_OH2 + norm_H12 * norm_H12 ) / ( 2 * norm_OH1 * norm_OH1 * norm_OH2 );
	double partial_cos_OH2 = ( - norm_OH1 * norm_OH1 + norm_OH2 * norm_OH2 + norm_H12 * norm_H12 ) / ( 2 * norm_OH1 * norm_OH2 * norm_OH2 );
	double partial_cos_H12 = ( - norm_H12 ) / ( norm_OH1 * norm_OH2 );

	double partial_u_OH1 = 2 * k_l * (norm_OH1 - standard_l) + 2 * k_theta * (theta - standard_theta) * partial_theta_cos * partial_cos_OH1;
	double partial_u_OH2 = 2 * k_l * (norm_OH2 - standard_l) + 2 * k_theta * (theta - standard_theta) * partial_theta_cos * partial_cos_OH2;
	double partial_u_H12 = 2 * k_theta * (theta - standard_theta) * partial_theta_cos * partial_cos_H12;

	VEC3 partial_OH1_vec = ( 1 / norm_OH1 ) * vec_OH1;
	VEC3 partial_OH2_vec = ( 1 / norm_OH2 ) * vec_OH2;
	VEC3 partial_H12_vec = ( 1 / norm_H12 ) * vec_H12;

	Force_O  = -1 * ( + partial_u_OH1 * partial_OH1_vec + partial_u_OH2 * partial_OH2_vec );
	Force_H1 = -1 * ( + partial_u_H12 * partial_H12_vec - partial_u_OH1 * partial_OH1_vec );
	Force_H2 = -1 * ( - partial_u_H12 * partial_H12_vec - partial_u_OH2 * partial_OH2_vec );
}

void Intramolecular_force(double natom, MAT3& coord, MAT3& force, VEC3& box_range, double& potential){
	VEC3 Oxyz, H1xyz, H2xyz;
	VEC3 Force_O, Force_H1, Force_H2;

	for (int index = 0; index < natom; index++){
		force[index] = {0, 0, 0};
	}

	double single_potential = 0;
	for (int index = 0; index < natom; index = index + 3){
		Oxyz  = coord[ index ];
		H1xyz = coord[ index + 1];
		H2xyz = coord[ index + 2];

		Intramolecular_force_for_single(Oxyz, H1xyz, H2xyz, box_range, Force_O, Force_H1, Force_H2, single_potential);
		force[ index ] = Force_O;
		force[ index + 1 ] = Force_H1;
		force[ index + 2 ] = Force_H2;

		potential += single_potential;
	}
}

void Load_Coord(MAT3& mat,string coord_path,int natom, int dim = 3){
    ifstream fxyz(coord_path);
    for (int i = 0; i < natom; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            fxyz >> mat[i][j];
        }
    }
    fxyz.close();
}

void load_init_file(double natom, MAT3& matrix, string loading_path){
	ifstream load_file(loading_path,ios::in);
	string line_buffer;
	int value_index = 0;
	VEC3 data;

	for(int index = 0;index <5;index++){
		load_file >> line_buffer;
	}

	while(load_file >> line_buffer){
		if(value_index % 4 == 0){
			value_index += 1;
			continue;
		}
		double single_data = atof(line_buffer.c_str());
		data.push_back(single_data);
		value_index += 1;
	}
	
	for (int index = 0; index < natom; index ++){
		for(int dim = 0; dim < 3; dim++){
			matrix[index][dim] = data[index * 3 + dim];
		}
	}
}

void VelocityInit(int natom, MAT3& Velocity)
{
	double KbT = kbp23 * nan23 * TargetTemp;
	mt19937 mt_rand(time(0));

    float vsigO = sqrt(KbT / ( UnitScale * MassO));
    float vsigH = sqrt(KbT / ( UnitScale *  MassH));
    normal_distribution<float> dO(0, vsigO);
    normal_distribution<float> dH(0, vsigH);
	for (int index = 0; index < natom; index += 3){
		for (int dim = 0; dim < 3; dim++){
			Velocity[index][dim] = dO(mt_rand);
		}
		for (int dim = 0; dim < 3; dim++){
			Velocity[index + 1][dim] = dH(mt_rand);
		}
		for (int dim = 0; dim < 3; dim++){
			Velocity[index + 2][dim] = dH(mt_rand);
		}
	}
}

void VelocityInit_Mol(int natom, MAT3& Velocity){
	double KbT = kbp23 * nan23 * TargetTemp;
	mt19937 mt_rand(time(0));
    float vsig = sqrt(KbT / ( UnitScale * (MassO + 2 * MassH)));
    normal_distribution<float> dMol(0, vsig);
	VEC3 velocity_buffer = {0,0,0};
	for (int index = 0; index < natom; index += 3){
		for (int dim = 0; dim < 3; dim++){
			velocity_buffer[dim] = dMol(mt_rand);
		}
		Velocity[index] = velocity_buffer;
		Velocity[index + 1] = velocity_buffer;
		Velocity[index + 2] = velocity_buffer; 
	}
}

void AcceleraInit(int natom, MAT3& Accelera){
	for( int index = 0; index < natom; index++){
		Accelera[index] = {0,0,0};
	}
} 

void ChargeInit(int natom, VEC3& Charge){
	for (int index = 0 ;index < nAtom; index = index + 3){
		Charge[index] = ChargeO;
		Charge[index + 1 ] = ChargeH;
		Charge[index + 2 ] = ChargeH;
	}
}

void MassInit(int natom, VEC3& Mass){
	for (int index = 0 ;index < nAtom; index = index + 3){
		Mass[index] = MassO;
		Mass[index + 1 ] = MassH;
		Mass[index + 2 ] = MassH;
	}
}

void Velocity_verlet(int natom, MAT3& position, MAT3& velocity, MAT3& accelera, int DeltaT){
	for (int index = 0; index < natom; index++){
		position[index] = position[index] + deltaT * velocity[index] + 0.5 * deltaT * deltaT * accelera[index];
		velocity[index] = velocity[index] + 0.5 * deltaT * accelera[index];
	}
}

void Velocity_verlet(int natom, MAT3& velocity, MAT3& accelera, int DeltaT){
	for (int index = 0; index < natom; index++){
		velocity[index] = velocity[index] + 0.5 * deltaT * accelera[index];
	}
}

void Calculate_force(MAT3& position, MAT3& Accelera, MAT3& force_inter, MAT3& force_intra,VEC3& Box, VEC3& Charge, VEC3& Mass, double& potential_inter, double& potential_intra, double lj_c0, double lj_c3, double lj_c4,double erfc_c0,double erfc_c3,double erfc_c4){
	potential_inter = 0;
	potential_intra = 0;
	Intermolecular_force(nAtom, Charge, position, force_inter, Box, potential_inter, b1_sh,b2_sh,sig_lj,eps_lj,alpha_eel,r4pie_eel,lj_c0,lj_c3,lj_c4,erfc_c0,erfc_c3,erfc_c4);
	Intramolecular_force(nAtom, position, force_intra, Box, potential_intra);
	for (int index = 0 ; index < nAtom; index++){
		Accelera[index] = (1 / Mass[index]) * (force_inter[index] + force_intra[index]);
	}
}

void Calculate_kinetic(int natom, MAT3& velocity, VEC3& Mass, double& kinetic, double& temp){
	kinetic = 0;
	double norm_buffer = 0;
	for( int index = 0; index < natom; index+= 1){
		norm_buffer = norm(velocity[index]);
		kinetic += 0.5 * Mass[index] * norm_buffer * norm_buffer;
	}
	temp = 2. * UnitScale * kinetic / (3. * natom * kbp23 * nan23);
}

void Flexible_MD(string init_coord_file, string save_info_folder, string init_velocity_file = "None"){

	MAT3 position(nAtom,VEC3(3));
	MAT3 velocity(nAtom,VEC3(3));
	MAT3 accelera(nAtom,VEC3(3));
	VEC3 Box = {Box_x, Box_y, Box_z};
	VEC3 Charge(nAtom);
	VEC3 Mass(nAtom);

	double temp_buffer;
	double potential, kinetic;
	double lj_c0,lj_c3,lj_c4;
	double eel_c0, eel_c3, eel_c4;

	// Initialize position
	// Load_Coord(position, init_coord_file, nAtom);
	load_init_file(nAtom, position, init_coord_file);
	// Initialize velocity
	VelocityInit(nAtom, velocity);
	// Initialize accelera
	AcceleraInit(nAtom, accelera);
	// Initialize charge
	ChargeInit(nAtom, Charge);
	// Initialize Mass
	MassInit(nAtom, Mass);
	// Initialize shift parameters
	shift_parameter_LJ(sig_lj,b1_sh,b2_sh,lj_c0,lj_c3,lj_c4);
	shift_parameter_erfc(alpha_eel,b1_sh,b2_sh,eel_c0,eel_c3,eel_c4);

	// Log basic parameters
	ofstream record_file(save_info_folder + "train_info.txt");
	time_t now;
	now = time(&now);
	char* curr_time = ctime(&now); 
	record_file << setw(25) << "# [ TIME START ] :\t" << curr_time;
	record_file << "# =====================    CONFIGURATION PARAMETERS    =====================" << endl;
	record_file << setw(25) << "# [ ATOM_NUMBER ] :\t" << nAtom << endl;
	record_file << setw(25) << "# [ BOX_LENGTH ] :\t" << Box_x << ", " << Box_y << ", " << Box_z << endl;
	record_file << "# =====================    ATOMIC PARAMETERS    =====================" << endl;
	record_file << setw(25) << "# [ ATOM_CHARGE ] :\t" << ChargeO << " , " << ChargeH << endl;
	record_file << setw(25) << "# [ ATOM_MASS ] :\t" << MassO << " , " << MassH << endl;
	record_file << "# =====================    PAIRWISE PARAMETERS    =====================" << endl;
	record_file << setw(25) << "# [ SHIFT_b1 ] :\t" << b1_sh << endl;
	record_file << setw(25) << "# [ SHIFT_b2 ] :\t" << b2_sh << endl;
	record_file << setw(25) << "# [ LJ_SIGMA ] :\t" << sig_lj << endl;
	record_file << setw(25) << "# [ LJ_EPSCLON ] :\t" << eps_lj << endl;
	record_file << setw(25) << "# [ LJ_SHIFT_C0 ] :\t" << lj_c0 << endl;
	record_file << setw(25) << "# [ LJ_SHIFT_C3 ] :\t" << lj_c3 << endl;
	record_file << setw(25) << "# [ LJ_SHIFT_C4 ] :\t" << lj_c4 << endl;
	record_file << setw(25) << "# [ R4PIE_EL ] :\t" << r4pie_eel << endl;
	record_file << setw(25) << "# [ ALPHA_EK ] :\t" << alpha_eel << endl;
	record_file << setw(25) << "# [ EEL_SHIFT_C0 ] :\t" << eel_c0 << endl;
	record_file << setw(25) << "# [ EEL_SHIFT_C3 ] :\t" << eel_c3 << endl;
	record_file << setw(25) << "# [ EEL_SHIFT_C4 ] :\t" << eel_c4 << endl;
	record_file << "# =====================    MD PARAMETERS    =====================" << endl;
	record_file << setw(25) << "# [ DELTA ] :\t" << deltaT << endl;
	record_file << setw(25) << "# [ TOTAL_STEP ] :\t" << Totalstep << endl;
	record_file << setw(25) << "# [ CHECK_STEP ] :\t" << mchk << endl;
	record_file << "# =====================    MD LOOP    =====================" << endl;
	record_file  << setw(20) << "# LOOP INDEX" << setw(20) << "TEMPRETURE" << setw(20) << "TOTAL ENERGY" << setw(20) << "KINETIC" << setw(20) << "POTENTIAL INTER" << setw(20)<< "POTENTIAL INTRA" << setw(20) << "VIRIAL" << endl;

	// Calculate force 
	MAT3 force_inter(nAtom,VEC3(3));
	MAT3 force_intra(nAtom,VEC3(3));

	double potential_inter = 0;
	double potential_intra = 0;

	Calculate_force(position,accelera,force_inter,force_intra,Box,Charge,Mass,potential_inter,potential_intra,lj_c0,lj_c3,lj_c4,eel_c0,eel_c3,eel_c4);
	Calculate_kinetic(nAtom, velocity, Mass, kinetic, temp_buffer);
	int step = 0;
	while (step < Totalstep)
	{
		Velocity_verlet(nAtom, position, velocity, accelera, deltaT);
		Calculate_force(position,accelera,force_inter,force_intra,Box,Charge,Mass,potential_inter,potential_intra,lj_c0,lj_c3,lj_c4,eel_c0,eel_c3,eel_c4);
		Velocity_verlet(nAtom, velocity, accelera, deltaT);
		if(step % mchk == 0){
			Calculate_kinetic(nAtom, velocity, Mass, kinetic, temp_buffer);
			record_file << setw(20) << step << setw(20) << temp_buffer << setw(20) << (kinetic + potential_intra + potential_inter)/nAtom << setw(20) << kinetic/nAtom << setw(20) << potential_inter/nAtom << setw(20) << potential_intra/nAtom << setw(20) << 0 << endl;
		}
		step += 1;
	}

}

int main(){
	string xyz_path_Hu = "/DATA/users/yanghe/projects/NeuralNetwork_PES/Identity_neural_network/Flexible_water_model/code/init.spce.xyz";
	string save_trj_path = "./";
	
	Flexible_MD(xyz_path_Hu,save_trj_path);
}