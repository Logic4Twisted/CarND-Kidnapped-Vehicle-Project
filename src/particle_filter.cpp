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
	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);
	for (int i = 0; i < num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;
	//cout << "yaw_rate = " << yaw_rate << endl;
	for (int i = 0; i < num_particles; i++) {
		Particle p = particles[i];
		//cout << "before " << p.x << ", " << p.x << ", " << p.theta << endl;
		double mean_x, mean_y, mean_theta;
		if (abs(yaw_rate) > 1.0e-10) {
			mean_x = velocity/yaw_rate*(sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
			mean_y = velocity/yaw_rate*(cos(p.theta) - cos(p.theta + yaw_rate*delta_t));
			mean_theta = yaw_rate*delta_t;
		}
		else {
			mean_x = velocity*delta_t*cos(p.theta);
			mean_y = velocity*delta_t*sin(p.theta);
			mean_theta = 0.0;
		}
		std::normal_distribution<double> dist_x(mean_x,std_pos[0]);
		std::normal_distribution<double> dist_y(mean_y,std_pos[1]);
		std::normal_distribution<double> dist_theta(mean_theta,std_pos[2]);
		Particle moved_particle;
		moved_particle.x = p.x + dist_x(gen);
		moved_particle.y = p.y + dist_y(gen);
		moved_particle.theta = normalize_angle(p.theta + dist_theta(gen));
		particles.erase(particles.begin());
		particles.push_back(moved_particle);
		//cout << "after " << moved_particle.x << ", " << moved_particle.x << ", " << moved_particle.theta << endl;
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	//cout << endl << "predicted: " << endl;
	//for (int j = 0; j < predicted.size(); j++) {
	//	cout << "> " << predicted[j].id << " (" << predicted[j].x << ", " << predicted[j].y << endl;
	//}
	cout << flush;
	for (int i = 0; i < observations.size(); i++) {
		double current_min = 2501.0;
		int index = -1;
		for (int j = 0; j < predicted.size(); j++) {
			if (index == -1 || dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y) < current_min) {
				current_min = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
				index = j;
			}
		}
		observations[i].id = index;
		//cout << "observed : (" << observations[i].x << ", " << observations[i].y << ") " << "\n";
		//cout << "matched : " << index << " (" << predicted[index].x << ", " << predicted[index].y << ") " << "\n";
		//cout << flush;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
	std::vector<LandmarkObs> landmarks;
	for (int i = 0; i < map_landmarks.landmark_list.size(); i++) {
		LandmarkObs l;
		l.x = map_landmarks.landmark_list[i].x_f;
		l.y = map_landmarks.landmark_list[i].y_f;
		l.id = map_landmarks.landmark_list[i].id_i;
		landmarks.push_back(l);
	}
	
	weights.clear();
	std::default_random_engine gen;
	for (int i = 0; i < particles.size(); i++) {
		double theta = particles[i].theta;
		std::vector<LandmarkObs> transf_observations;
		for (int j = 0; j < observations.size(); j++) {
			double x = observations[j].x;
			double y = observations[j].y;
			LandmarkObs l;
			l.x = particles[i].x + x*cos(theta) - y*sin(theta);
			l.y = particles[i].y + x*sin(theta) + y*cos(theta);
			transf_observations.push_back(l);
		}
		dataAssociation(landmarks, transf_observations);
		particles[i].associations.clear();
		particles[i].sense_x.clear();
		particles[i].sense_y.clear();
		double weight = 1.0;
		for (int j = 0; j < transf_observations.size(); j++) {
			int index = transf_observations[j].id;
			//cout << "index = " << index << "\n";
			//cout << flush;
			//cout << "(" << transf_observations[j].x << ", " << transf_observations[j].y << ") - (" << landmarks[index].x << ", " << landmarks[index].y << ")" << "\n";
			weight *= calc_prob(transf_observations[j].x, landmarks[index].x, std_landmark[0]);
			weight *= calc_prob(transf_observations[j].y, landmarks[index].y, std_landmark[1]);
			//cout << " - " << weight << " - \n";
			//particles[i].associations.push_back(index);
			//particles[i].sense_x.push_back(transf_observations[j].x);
			//particles[i].sense_y.push_back(transf_observations[j].y);
		}
		//cout << i << " weight =" << weight << "\n";
		weights.push_back(weight);
		particles[i].weight = weight;
	}
	//cout << "---------------------------------";
	//cout << flush;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::default_random_engine gen;
	std::discrete_distribution<int> dist(weights.begin(),weights.end());
	std::vector<Particle> resampled;
	for (int i = 0; i < num_particles; i++) {
		resampled.push_back(particles[dist(gen)]);
	}
	particles = resampled;
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

void ParticleFilter::printParticles() {
	cout << "Paricles :" << endl;
	for (int i = 0; i < num_particles; i++) {
		Particle p = particles[i];
		std:cout << "(" << p.x << ", " << p.y << ", " << p.theta << ")" << std::endl;
	}
	cout << endl;
}
