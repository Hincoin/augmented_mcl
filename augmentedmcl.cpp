// AugmentedMCL.cpp : Defines the entry point for the console application.
//
// Written by: Alejandro Lucena
// Re-localization of robot after kidnapping
// Particle Filter


#include <iostream>
#include <algorithm>
#include <ctime>
#include <unordered_set>
#include <random>
#include <vector>
#include <algorithm>
#include <cmath>
#define PI 3.14159265358979323846
#define TWO_PI  2*PI
#define ALPHA_SLOW  .05
#define ALPHA_FAST .5
using namespace std;
vector<vector<double>> landmarks;
class particle
{
public:	
	double x, y,theta,range_noise, bearing_noise ,forward_noise, rotation_noise;
	double importance_weight;
	vector<pair<double,double>> measurements;
	particle(double x1, double y1, double t, double range_n, double bn, double fn, double rn)
	{
		x = x1; y = y1; theta = t;
		range_noise = range_n; bearing_noise = bn; forward_noise = fn; rotation_noise = rn;
		importance_weight = 1.0;
	}

	friend ostream& operator << (ostream& os,const particle& p)
	{
		os << "x: " << p.x << " y: " << p.y << " theta: " << p.theta;
		return os;
	}

	
	
};

// Global Variables
double range_sig= 1;
double for_sig = .2;
double rot_sig = .2;
double bear_sig = 1;
particle robot(7.0,7.0,PI/4,range_sig,bear_sig,for_sig,rot_sig);


inline vector<double> get_data(const particle& ref, const vector<double>& lmark, int noise)
{
	std::default_random_engine generator; 
	std::normal_distribution<double> range_dist(0.0,ref.range_noise);
	std::normal_distribution<double> bear_dist(0.0,ref.bearing_noise);
	vector<double> data;
	data.push_back(sqrt(pow(ref.x - lmark[0],2) + pow(ref.y - lmark[1],2))+ noise*  range_dist(generator));
	data.push_back(atan2(lmark[1]-ref.y,lmark[0]-ref.x) + noise* bear_dist(generator));
	return data;
}

vector<particle> make_particles(int M,double range_sig,double bearing_sig,double for_sig, double rot_sig)
{
	std::default_random_engine generator; 
	/*std::normal_distribution<double> sens_dist(0.0,sens_sig);
	std::normal_distribution<double> forward_dist(0.0,for_sig);
	std::normal_distribution<double> rot_dist(0.0,rot_sig);
	*/
	std::uniform_real_distribution<double> position_dist(0,100.0);
	std::uniform_real_distribution<double> theta_dist(0,2 * PI);
	vector<particle> particles;
	for(int i =0; i < M; ++i)
	{
		//particle p(position_dist(generator),position_dist(generator),theta_dist(generator),sens_dist(generator),forward_dist(generator),rot_dist(generator));
		particle p(position_dist(generator),position_dist(generator),theta_dist(generator),range_sig,bearing_sig,for_sig,rot_sig);
		particles.push_back(p);
	}
	return particles;

}


void motion(particle& p, double distance, double turn_angle)
{


	std::default_random_engine generator; 
	std::normal_distribution<double> forward_dist(0.0,p.forward_noise);
	std::normal_distribution<double> rot_dist(0.0,p.rotation_noise);
	turn_angle += rot_dist(generator);
	distance += forward_dist(generator);

	double beta = (distance / 20.0) * tan(turn_angle);
	if(abs(beta) <= .001)
	{
		p.x += distance * cos(p.theta);
		p.y += distance * sin(p.theta);
		p.theta += beta;
		p.theta = fmod(p.theta,TWO_PI);

	}
	else
	{
		double R = distance / beta;
		double cx = p.x - R * sin(p.theta);
		double cy = p.y + R* cos(p.theta);
		p.theta += beta;
		p.theta = fmod(p.theta,TWO_PI);
		p.x = cx + R * sin(p.theta);
		p.y = cy -  R* cos(p.theta);


	}
}

vector<double> make_vector(double t1, double t2)
{
	vector<double> data;
	data.push_back(t1);
	data.push_back(t2);
	return data;
}
inline void importance_weights(vector<particle>& particles, const vector<pair<double,double>>& real_data, double& slow, double& fast)
{
	double total_sum = 0.0;
	double avg = 0;
	int M = particles.size();
	for(auto& p : particles)
	{
		for(const auto& landmark : landmarks)
		{
			auto data_vec = get_data(p,landmark,0);
			p.measurements.push_back(make_pair(data_vec[0],data_vec[1]));
		}
		int i =0 ;
		for(const auto& pr : p.measurements)
		{

			double expon_range = exp(-1 *  pow(pr.first-real_data[i].first,2) / (2 * p.range_noise));
			expon_range /= (sqrt(2*PI*p.range_noise));
			double expon_bear = exp(-1 * pow(pr.second-real_data[i].second,2) / (2  * p.bearing_noise));
			expon_bear /= (sqrt(2 * PI  * p.bearing_noise));
			p.importance_weight *= (expon_range + expon_bear);
			++i;
			
			
		}
		total_sum += p.importance_weight;
		avg += (p.importance_weight/M);
	}

	slow += ALPHA_SLOW * (avg-slow);
	fast += ALPHA_FAST * (avg-fast);
	for(auto& p  :particles)
	{
		p.importance_weight /= total_sum;
	}
	
}
vector<particle> resample(const vector<particle>& prior,double slow, double fast,int correspondence)
{
	auto r = robot;
	vector<particle> posterior;
	vector<double> importance_weights;
	for(const auto& x : prior)
		importance_weights.push_back(x.importance_weight);
	std::uniform_real_distribution<double> add_random(0,1);
	std::uniform_real_distribution<double> yeta_dist(0,2*PI);
	std::normal_distribution<double> range_dist(0.0,prior[0].range_noise);
	std::normal_distribution<double> bear_dist(0.0,prior[0].bearing_noise);
	std::default_random_engine generator;
	double vl = 0;
	double threshold = 1.0 - (fast/slow);
	cout << "\nthreshold = " << threshold << "\n\n";
	if(threshold > vl) vl = threshold;
	int index = rand() % prior.size();
	double beta =0.0;

	double max_val =*max_element(importance_weights.begin(),importance_weights.end()); 

	for(int i =0 ; i < prior.size(); ++i)
	{
		if(add_random(generator) < vl)
		{
			// add random pose according to the circle of possible poses around the landmark
			
			double yeta = yeta_dist(generator);
			double range_prime = r.measurements[correspondence].first + range_dist(generator);
			double bear_prime = r.measurements[correspondence].second + bear_dist(generator);
			double x_prime = landmarks[correspondence][0] + (range_prime * cos(yeta));
			double y_prime = landmarks[correspondence][1] + (range_prime * sin(yeta));
			double theta_prime = yeta - PI - bear_prime;
			if(theta_prime < 0)
				theta_prime += 2*PI;
			particle sampled(x_prime,y_prime,theta_prime,prior[i].range_noise,prior[i].bearing_noise,prior[i].forward_noise,prior[i].rotation_noise);
			posterior.push_back(sampled);
		}
		else
		{
			beta += fmod(rand() ,( 2 * max_val));
			while(importance_weights[index] < beta)
			{
				beta -= importance_weights[index];
				++index;
				index %= prior.size();
			}
			particle new_particle = prior[index];
			new_particle.importance_weight=1.0;
			new_particle.measurements = vector<pair<double,double>>();
			posterior.push_back(new_particle);
			

		}
	}
	return posterior;
	
}
pair<double,double> average_particles(const vector<particle>& particles)
{
	double x_total = 0.0, y_total = 0.0;
	for(const auto& p : particles)
	{
		x_total += p.x;

		y_total += p.y;
	}
	x_total /= particles.size();
	y_total /= particles.size();
	return make_pair(x_total,y_total);
}
int main()
{
	srand(time(NULL));
	landmarks.push_back(make_vector(5.0,5.0));
	landmarks.push_back(make_vector(10.0,20.0));
	landmarks.push_back(make_vector(32.0,16.0));
	landmarks.push_back(make_vector(20.0,17.0));
	landmarks.push_back(make_vector(40.0,40.0));
	landmarks.push_back(make_vector(23.0,10.0));
	landmarks.push_back(make_vector(8.0,22.0));
	landmarks.push_back(make_vector(16.0,17.0));
	 
	int M = 1000;
	
	auto particles = make_particles(M,range_sig,bear_sig,for_sig,rot_sig);
	double slow = 0.0, fast = 0.0;
	for(int i = 0; i < 30; ++i)
	{
		motion(robot,sqrt(2),0.0);
		for(auto& p : particles)
			motion(p,sqrt(2),0.0);
		int correspond =landmarks.size()-1; 
		for(int x= 0; x < landmarks.size(); ++x)
		{
			auto vec = get_data(robot,landmarks[x],0);
			double dst = vec[0];//sqrt(pow(robot.x - landmarks[x][0],2) + pow(robot.y-landmarks[x][1],2));
			double bear = vec[1];//atan2(robot.y-landmarks[x][1],robot.x-landmarks[x][0]);
			robot.measurements.push_back(make_pair(dst,bear)); // robot perceives the environment
		}
		
		importance_weights(particles,robot.measurements,slow,fast);
		auto new_particles = resample(particles,slow,fast,correspond);
		particles = new_particles;
		cout << "robot = " << robot << "\n";
		auto guess = average_particles(particles);
		cout << "guess: " << guess.first << "    " << guess.second << "\n";
		int q= 10;
		if(i == 5)
		{
			robot.x = 50;
			robot.y = 50;
		}
		robot.measurements = vector<pair<double,double>>();

	}
	

	int br;

	cin >> br;
	return 0;
}



