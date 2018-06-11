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
	default_random_engine gen;
	num_particles = 100;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	particles.reserve(num_particles);
	weights.reserve(num_particles);
	
	for(int i=0;i<num_particles;i++){
	  particles[i].id = i;
	  particles[i].x = dist_x(gen);
	  particles[i].y = dist_y(gen);
	  particles[i].theta = dist_theta(gen);
	  particles[i].weight = 1.;
	  weights[i] = 1.;
	}
	  
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	//avoid division by zero
	default_random_engine gen;
	
	for(int i=0;i<num_particles;i++){
	
	  if (fabs(yaw_rate) > 0.001) {
	    particles[i].x += velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
		particles[i].y += velocity/yaw_rate*(-cos(particles[i].theta+yaw_rate*delta_t)+cos(particles[i].theta));
		particles[i].theta += yaw_rate*delta_t;
	  } else {
	    particles[i].x += velocity*delta_t*cos(particles[i].theta);
	    particles[i].y += velocity*delta_t*sin(particles[i].theta);
	  }
	  normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
	  normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
	  normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
	  
	  particles[i].x = dist_x(gen);
	  particles[i].y = dist_y(gen);
	  particles[i].theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	/* NOT USED
	
	double distance;
	double tmp_dist;
	for(unsigned int j=0,j<observations.size(),j++){
	  distance = 1e6;
	  for(unsigned int i=0;i<predicted.size();i++){
		tmp_dist = sqrt(pow(predicted(i).x-observations[j].x,2)+pow(predicted(i).y-observations[j].y,2));
		if(tmp_dist<distance){
		  distance = tmp_dist;
		  observations[j].id = i;
		}
	  }
	}
	
	*/
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
	
	//initialize variables
	const double sx = std_landmark[0];
	const double sy = std_landmark[1];
	
	const double xyden = 2*M_PI*sx*sy;
	const double xden = 2*sx*sx;
	const double yden = 2*sy*sy;
	
	double Gauss_dist;
	double obs_map_x, obs_map_y;
	double particle_distance;
	
	int argmin;
	
	vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
    vector<double> landmark_distances(landmarks.size());
	
	//variable initialization end
	
	for(int i=0;i<num_particles;i++){
	  
	  //particles[i].associations.clear();
	  //particles[i].sense_x.clear();
	  //particles[i].sense_y.clear();
	  
	  Gauss_dist = 1.;

	  for(unsigned int j=0;j<observations.size();j++){
	    
		// transform observations
        obs_map_x = observations[j].x*cos(particles[i].theta)-observations[j].y*sin(particles[i].theta)+particles[i].x;
        obs_map_y = observations[j].x*sin(particles[i].theta)+observations[j].y*cos(particles[i].theta)+particles[i].y;
		
		//particles[i].sense_x.push_back(obs_map_x);
		//particles[i].sense_y.push_back(obs_map_y);
		
		fill(landmark_distances.begin(),landmark_distances.end(),sensor_range*10.);
	  
		for (unsigned int k=0;k<landmarks.size();k++) {
        
		  // within sensor range
          //particle_distance = sqrt(pow(particles[i].x-landmarks[k].x_f,2)+pow(particles[i].y-landmarks[k].y_f,2));
          //if (particle_distance <= sensor_range)
            landmark_distances[k] = sqrt(pow(obs_map_x-landmarks[k].x_f,2)+pow(obs_map_y-landmarks[k].y_f,2));
	    }
		// nearest neighbor
        argmin = distance(landmark_distances.begin(),min_element(landmark_distances.begin(),landmark_distances.end()));
		//particles[i].associations.push_back(argmin);
      
        // Multi-variate Gaussian distribution
        Gauss_dist *= exp(-pow(obs_map_x-landmarks[argmin].x_f,2)/xden-pow(obs_map_y-landmarks[argmin].y_f,2)/yden)/xyden;
	  }
	  particles[i].weight = Gauss_dist;
	  weights[i] = Gauss_dist;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
	// Vector for resampled particles
    vector<Particle> resampled_particles(num_particles);
  
    // discrete distribution resamples particles by weight
    default_random_engine gen;
  
    for (int i=0;i<num_particles;i++) {
      discrete_distribution<int> idx(weights.begin(), weights.end());
      resampled_particles[i] = particles[idx(gen)];
    }
	
    // Resample
    particles = resampled_particles;
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
