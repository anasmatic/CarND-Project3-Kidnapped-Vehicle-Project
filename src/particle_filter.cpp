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
	cout << "init" << endl;
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 10;//random number by trial and error
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		double sample_x, sample_y, sample_theta;
		// Sample  and from these normal distrubtions like this: 
		//	sample_x = dist_x(gen);
		//	where "gen" is the random engine initialized earlier.
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1;
		particle.id = i;

		particles.push_back(particle);
		weights.push_back(1);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	cout << "prediction" << endl;
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//θ​f​​ = θ​0​​ + ​θ​˙​​(dt)
	default_random_engine gen;
	double velocity_over_yaw_rate = (yaw_rate == 0) ? velocity / 1e-6 : velocity/yaw_rate ;
	
	//Particle particle; for (int i = 0; i < particles.size(); i++) { particle = particles[i]
	for (Particle& particle : particles) {
		double x = 0, y = 0 , theta = 0;
		double yaw_rate_time_delta_t = yaw_rate*delta_t;
		x = particle.x + (velocity_over_yaw_rate * ( sin(particle.theta+(yaw_rate_time_delta_t)) - sin(particle.theta)));//xf = x​0​​ + v/​​θ​˙​[sin(θ​0​​ + ​θ​˙​​(dt))−sin(θ​0​​)]
		y = particle.y + (velocity_over_yaw_rate * ( cos(particle.theta) - cos(particle.theta+(yaw_rate_time_delta_t))));//y​f​​ = y​0​​ + v/​​θ​˙​[cos(θ​0​​)−cos(θ​0​​ + ​θ​˙​​(dt))]
		theta = particle.theta + (yaw_rate * delta_t);
		
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
	cout << "dataAssociation" << predicted.size() <<","<< observations.size() << endl;
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	/*for (LandmarkObs obs : observations)
	{
//cout << "		predicted.size " << predicted.size() << endl;
		double min_dist = 9999999999;
		for (LandmarkObs pred : predicted)
		{
			double distance = sqrt(pow(pred.x - obs.x, 2) + pow(pred.y - obs.y, 2));
//cout << "		distance " << distance << endl;
			if (distance <= min_dist) {
//cout << "			obs.id " << obs.id << ", pred.id "<< pred.id << endl;
				obs.id = pred.id;
			}
		}
	}*/
	for (int i = 0; i < observations.size(); ++i) {
		double min_distance = 1000;
		int closest_landmark = -1;
		for (int j = 0; j < predicted.size(); ++j) {
			double distance = sqrt((predicted[j].x - observations[i].x)*(predicted[j].x - observations[i].x) +
				(predicted[j].y - observations[i].y)*(predicted[j].y - observations[i].y));
			if (distance < min_distance) {
				min_distance = distance;
				closest_landmark = j;
			}
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	cout << "updateWeights "<< particles.size() << endl;
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	double weights_sum = 0;
	//for (Particle& particle : particles) {
	Particle particle; for (int i = 0; i < particles.size(); i++) {
		particle = particles[i];
		// transform Observations -in car coordinate system- into "map coordinates"
		vector<LandmarkObs> observations_map_coor;
		for (LandmarkObs const& obs : observations) {
			LandmarkObs obs_map;
			obs_map.x = particle.x + (cos(particle.theta) * obs.x) - (sin(particle.theta) * obs.y);
			obs_map.y = particle.y + (sin(particle.theta) * obs.x) + (cos(particle.theta) * obs.y);
			obs_map.id = particle.id;
			observations_map_coor.push_back(obs_map);
		}
		//add land mark if within sensor range
		vector<LandmarkObs> landmarks_inrange;
		for (Map::single_landmark_s const& landmark : map_landmarks.landmark_list)
		{
			double distance = sqrt(pow(landmark.x_f - particle.x, 2) +pow(landmark.y_f - particle.y, 2));
			if (distance <= sensor_range) {
				LandmarkObs obs;
				obs.x = landmark.x_f; 
				obs.y = landmark.y_f;
				obs.id = landmark.id_i;
				landmarks_inrange.push_back(obs);
			}
		}

		dataAssociation(landmarks_inrange, observations_map_coor);
		particle.weight = 1;
		cout << "	- "<< particle.weight << endl;
		for (LandmarkObs obs : observations_map_coor) {
			LandmarkObs map_landmark = landmarks_inrange[obs.id];
			double std_x = std_landmark[0], std_y = std_landmark[1];
			cout << "		- std x:" << std_x <<" std y:" << std_y << endl;
			double step1 = 1 / (2 * M_PI*std_x*std_y);
			cout << "		- step1:" << step1 << endl;
			double e_pow_step2_1 = pow((obs.x - map_landmark.x),2) / (2 * pow(std_x,2));
			double e_pow_step2_2 = pow((obs.y - map_landmark.y),2) / (2 * pow(std_y,2));
			double e_pow = (e_pow_step2_1 + e_pow_step2_2) * -1;
			cout << "		- e_pow: " << e_pow_step2_1 <<"+"<< e_pow_step2_2<<" *-1 = "<< e_pow << endl;
			cout << "		- res			:" << step1 * exp(e_pow) << endl;
			particle.weight *= (step1 * exp(e_pow));
			cout << "			-*particle.weight:" << particle.weight << endl;
		}
		cout << "	-- "<< particle.weight << endl;
		weights[i] = particle.weight;
		weights_sum += particle.weight;
	}
	cout << "	end loop. sum:"<< weights_sum << endl;
	//normlize weights
	if (weights_sum != 0) 
		for (int i = 0; i < particles.size(); i++) {
			particles[i].weight = particles[i].weight / weights_sum;
			weights[i] = particles[i].weight;
		}
	cout << " end norm" << endl;
}

void ParticleFilter::resample() {
	cout << "resample" << endl;
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	double max_weight = -999999;
	for (Particle& particle : particles)
		if (particle.weight > max_weight)
			max_weight = particle.weight;
	double mw = max_weight * 2;
	vector<Particle> resample_particles;
	//because I don't know how to utilize discrete_distribution with resample wheel
	std::default_random_engine generator;
	std::uniform_int_distribution<int> dist_num_particles(0, num_particles);
	int index = dist_num_particles(generator);
	std::uniform_int_distribution<int> dist_mw(0, mw);
	double beta = 0.0;
	for (int i = 0; i < num_particles; i++)
	{
		beta += dist_mw(generator);
		while (weights[index] < beta) {
			beta -= weights[index];
			index = (index + 1) % num_particles; // %N to convert 1000 to 0
		}
		resample_particles.push_back(particles[index]);
	}
	particles = resample_particles;
}

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
