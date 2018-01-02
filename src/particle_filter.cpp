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
default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 40;

	// This creates a normal (Gaussian) distribution for x,y and theta.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// generating particles
	for (int i = 0; i < num_particles; i++) {
		Particle sample_particle;

		sample_particle.id = i;
		sample_particle.x = dist_x(gen);
		sample_particle.y = dist_y(gen);
		sample_particle.theta = dist_theta(gen);
		sample_particle.weight = 1.0;

		particles.push_back(sample_particle);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Predicting particle location using motion model
	for (int i = 0; i < particles.size(); i++) {
		default_random_engine gen;

		if (yaw_rate>0.0001){
			// predict x, y and theta
			particles[i].x = particles[i].x + (velocity / yaw_rate)*(sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + (velocity / yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			particles[i].theta = particles[i].theta + yaw_rate * delta_t;
		} else {
			// predict x, y and theta
			particles[i].x = particles[i].x + velocity * cos(particles[i].theta) * delta_t;
			particles[i].y = particles[i].y + velocity * sin(particles[i].theta) * delta_t;

		}

		// add noise
		// x
		normal_distribution<double> pred_noise_x(particles[i].x, std_pos[0]);
		particles[i].x = pred_noise_x(gen);

		normal_distribution<double> pred_noise_y(particles[i].y, std_pos[1]);
		particles[i].y = pred_noise_y(gen);

		normal_distribution<double> pred_noise_theta(particles[i].theta, std_pos[2]);
		particles[i].theta = pred_noise_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	double calc_diff;

	for(int i=0; i<observations.size();i++){
		double min_diff = 1.0e99;
		int aux_id = -1;

		for(int j=0; j<predicted.size(); j++){
			calc_diff = sqrt(pow(predicted[j].x - observations[i].x,2) + pow(predicted[j].y - observations[i].y,2));

			if(calc_diff<min_diff){
				min_diff = calc_diff;
				aux_id = predicted[j].id;
			}
		}
		observations[i].id = aux_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	// stepping through each particle
	for (int i = 0; i < particles.size(); i++) {
		vector<LandmarkObs> obs;
		vector<LandmarkObs> near_landm;
		obs = observations;
		double x_m;
		double y_m;

		// transforming coordinates of observations(measurements) of car to each particle
		for (int j = 0; j < obs.size(); j++){
			x_m = particles[i].x + cos(particles[i].theta) * obs[j].x - sin(particles[i].theta) * obs[j].y;
			y_m = particles[i].y + sin(particles[i].theta) * obs[j].x + cos(particles[i].theta) * obs[j].y;
			obs[j].x = x_m;
			obs[j].y = y_m;
		}

		// extracting landmarks within sensor range
		for (int k = 0; k < map_landmarks.landmark_list.size(); k++){
			double d;

			d = sqrt(pow((map_landmarks.landmark_list[k].x_f - particles[i].x),2) + pow((map_landmarks.landmark_list[k].y_f - particles[i].y),2));

			if (d <= sensor_range){
				near_landm.push_back(LandmarkObs{map_landmarks.landmark_list[k].id_i ,map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f});
			}
		}

		// Associating landmarks to observations
		dataAssociation(near_landm, obs);

		// stepping through observations to calculate their weight and eventually multiply all observation weights to
		// get each particle weight
		for (int j = 0; j < obs.size(); j++){
			double o_x, o_y, s_x, s_y, mu_x, mu_y, obs_weight;
			o_x = obs[j].x;
			o_y = obs[j].y;
			s_x = std_landmark[0];
			s_y = std_landmark[1];

			for (int k = 0; k < near_landm.size(); k++){
				if (obs[j].id = near_landm[k].id){
					mu_x = near_landm[k].x;
					mu_y = near_landm[k].y;
				}
			}
			// calculate weight using normalization terms and exponent
			obs_weight =( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(o_x-mu_x,2)/(2*pow(s_x, 2)) + (pow(o_y-mu_y,2)/(2*pow(s_y,2))) ) );
			particles[i].weight *= obs_weight;
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> resam;
	double beta = 0.0;
	vector<double> extract_weights;
	double mw;

  // generate random starting index for resampling wheel
  uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

	// extract weights from particles
	for (int i = 0; i < particles.size(); i++){
		extract_weights.push_back(particles[i].weight);
	}

	// extracting max weight
	mw = *max_element(extract_weights.begin(), extract_weights.end());


  // uniform random distribution [0.0, max_weight)
  uniform_real_distribution<double> unirealdist(0.0, mw);

	for (int i = 0; i < particles.size(); i++){
		beta += unirealdist(gen) * 2.0;
		while(beta > extract_weights[index]){
			beta -= extract_weights[index];
			index = (index + 1) % particles.size();
		}
		resam.push_back(particles[index]);
	}
	particles = resam;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
