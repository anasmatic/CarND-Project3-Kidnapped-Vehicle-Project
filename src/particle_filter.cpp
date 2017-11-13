/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 75;//random number by trial and error
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		// Sample  and from these normal distrubtions like this: 
		//	sample_x = dist_x(gen);
		//	where "gen" is the random engine initialized earlier.
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1;
		particle.id = i;////
		particles.push_back(particle);
		weights.push_back(1);////
	}
	is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	//######## REQUIRED 1 #####################
	//update yaw_rate : when yaw _ate is ~equal to zero. we devide velocity over 0.001, else devide velocity over yaw_rate
	double velocity_over_yaw_rate = (fabs(yaw_rate) < 0.00001) ? velocity / 1e-3 : velocity / yaw_rate;
	//so instead of writing the same equation twice inside the foor loop,
	//	 once to handle zero division, and the other to handle normal situation
	//	I noticed that the only operation will be effected is (velocity/yaw_rate)
	//	so I calculated it here, just before we enter the loop
	//######## ######### #####################
	for (Particle& particle : particles) {
		//I don't need to handle division by zero in the for loop because it is already handled before it : velocity_over_yaw_rate
		double x = 0, y = 0, theta = 0;
		//do some calculations first to make the equations more readable
		double yaw_rate_time_delta_t = yaw_rate*delta_t;
		double theta_plus_yaw_rate_times_delta_t = particle.theta + (yaw_rate_time_delta_t);
		//calculate x, y ,theta
		x = particle.x + (velocity_over_yaw_rate * (sin(theta_plus_yaw_rate_times_delta_t) - sin(particle.theta)));//xf = x​0​​ + v/​​θ​˙​[sin(θ​0​​ + ​θ​˙​​(dt))−sin(θ​0​​)]
		y = particle.y + (velocity_over_yaw_rate * (cos(particle.theta) - cos(theta_plus_yaw_rate_times_delta_t)));//y​f​​ = y​0​​ + v/​​θ​˙​[cos(θ​0​​)−cos(θ​0​​ + ​θ​˙​​(dt))]
		theta = particle.theta + (yaw_rate * delta_t);	//θ​f​​ = θ​0​​ + ​θ​˙​​(dt)
		//add gaussian noise , with mean = updated partical pos values , and standerd deviation = std_pos[]
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[2]);
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i < observations.size(); ++i)
	{
		double min_dist = 9999999999;
		//for (LandmarkObs pred : predicted)
		for (int j = 0; j < predicted.size(); j++)
		{
			double distance = sqrt(pow(predicted[j].x - observations[i].x, 2) + pow(predicted[j].y - observations[i].y, 2));
			if (distance < min_dist) {
				min_dist = distance;
				observations[i].id = j;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
	const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
		//for (Particle& particle : particles) { 
	for (int i = 0; i < particles.size(); i++) {
		// transform Observations -in car coordinate system- into "map coordinates"
		vector<LandmarkObs> observations_map_coor;
		for (LandmarkObs const& obs : observations) {
			LandmarkObs obs_map;
			obs_map.x = particles[i].x + (cos(particles[i].theta) * obs.x) - (sin(particles[i].theta) * obs.y);
			obs_map.y = particles[i].y + (sin(particles[i].theta) * obs.x) + (cos(particles[i].theta) * obs.y);
			obs_map.id = particles[i].id;
			observations_map_coor.push_back(obs_map);
		}
		//add landmark if within sensor range
		vector<LandmarkObs> landmarks_inrange;
		for (Map::single_landmark_s const& landmark : map_landmarks.landmark_list)
		{
			double distance = sqrt(pow(landmark.x_f - particles[i].x, 2) + pow(landmark.y_f - particles[i].y, 2));
			if (distance <= sensor_range) {
				LandmarkObs obs;
				obs.x = landmark.x_f;
				obs.y = landmark.y_f;
				obs.id = landmark.id_i;
				landmarks_inrange.push_back(obs);
			}
		}

		dataAssociation(landmarks_inrange, observations_map_coor);

		//update weights according to how much partical match observation landmarks
		particles[i].weight = 1;
		for (LandmarkObs obs : observations_map_coor) {
			LandmarkObs map_landmark = landmarks_inrange[obs.id];
			double std_x = std_landmark[0], std_y = std_landmark[1];
			double step1 = 1 / (2 * M_PI*std_x*std_y);
			double e_pow_step2_1 = pow((obs.x - map_landmark.x), 2) / (2 * pow(std_x, 2));
			double e_pow_step2_2 = pow((obs.y - map_landmark.y), 2) / (2 * pow(std_y, 2));
			double e_pow = (e_pow_step2_1 + e_pow_step2_2) * -1;
			particles[i].weight *= (step1 * exp(e_pow));
		}
		weights[i] = particles[i].weight;
	}
}

//######## REQUIRED 2 #####################
void ParticleFilter::resample() {
	vector<Particle> resample_particles;
	// get all of the current weights
	vector<double> weights;
	for (int i = 0; i < num_particles; i++) {
	  weights.push_back(particles[i].weight);
	}
	// generate random starting index for resampling wheel
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());//used for index
	// get max weight
	double max_weight = *max_element(weights.begin(), weights.end());
	// uniform random distribution [0.0, max_weight)
	uniform_real_distribution<double> unirealdist(0.0, max_weight);//used for beta

	//assign to index rand value from weights.begin() to weights.end()
	int index = distribution(gen);
	double beta = 0.0;
	// spin the resample wheel!
	for (int i = 0; i < num_particles; i++) {
		beta += unirealdist(gen) * 2.0;//assign rand value from 0 to max_weight * 2
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles; // %N to convert 1000 to 0
		}
		resample_particles.push_back(particles[index]);
	}
	
	particles = resample_particles;
	//thanks you, resampling wheel is working now <3
}
//######## ######### #####################

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
