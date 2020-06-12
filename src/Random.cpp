#include <algorithm>
#include <cassert>
#include "Random.hpp"

std::random_device RandomDist::rd;
std::mt19937 RandomDist::rng = std::mt19937(RandomDist::rd());


RandomDist::RandomDist(double m, double s, int ns, bool n) 
  : mean(m), sd(s), nsample(ns), normdist(n)
{
    if (s <= 0) {
        throw(1);
    }
    if (ns <= 0) {
        throw(2);
    }
}


std::vector<double> RandomDist::generate_numbers() {
    std::vector<double> result(nsample);
    if (normdist) {
        normal(result);
    }
    else {
        uniform(result);
    }
    return result;
}


int RandomDist::binomial(int n, double p) {
	std::binomial_distribution<int> dbinom(n, p);
	return dbinom(rng);
}


int RandomDist::uniformIntSingle(int min, int max) {
	// init random distribution
	std::uniform_int_distribution<int> distr(min, max);
	
	// get one value
	return distr(RandomDist::rng);
}


double RandomDist::uniformDoubleSingle(double min, double max) {
	// init random distribution
	std::uniform_real_distribution<> distr(min, max);
	
	// get one value
	return distr(RandomDist::rng);
}


void RandomDist::uniformIntVector(std::vector<int>& toFill, int min, int max) {
	// init random distribution
	std::uniform_int_distribution<int> distr(min, max);

	// fill table with generated values
	std::generate(
		toFill.begin(),
		toFill.end(), 
		[&]() { 				
			return distr(RandomDist::rng); 
		}
	); 
}


void RandomDist::uniformDoubleVector(std::vector<double>& toFill, double min, double max) {
	// init random distribution
	std::uniform_real_distribution<> distr(min, max);

	// fill table with generated values
	std::generate(
		toFill.begin(),
		toFill.end(), 
		[&]() { 				
			return distr(RandomDist::rng); 
		}
	); 
}


void RandomDist::uniform(std::vector< double > &res) {
    double delta = sd * sqrt(3.0);
    double lower = mean - delta, upper = mean + delta;
    std::uniform_real_distribution<> unif(lower, upper);
    
    for (auto I = res.begin(); I != res.end(); I++) {
        *I = unif(rng);
    }
}


void RandomDist::normal(std::vector< double > &res) {
    std::normal_distribution<> norm(mean, sd);
   
    for (auto I = res.begin(); I != res.end(); I++) {
        *I = norm(rng);
    }
}


void RandomDist::multinomial(std::vector<unsigned int>& pop) {
	int n = 0;
	for (auto& count : pop)
		n += count;
		
	RandomDist::multinomial(pop, n);
}


void RandomDist::multinomial(std::vector<unsigned int>& pop, int n) {
	int total = 0;
	for (auto& count : pop)
		total += count;
	
	for (auto& count : pop) {
		
		// remaining parent population size should be 0 or more	
		assert(total >= 0);
		if (total == 0) {
			// the only way for the parent population to be 0 is if the 
			// current allele is not present anymore (wiped out)
			assert(count == 0);
			assert(n == 0);
			
			// if the parent population is empty, there are no more 
			// alleles to be distributed, thus we can exit the loop
			break;
		} 
		
		// generate new allele copy number
		double p = count * 1.0 / total;
		
		// reduce residual "gene pool"
		total -= count;
		
		// generate new number of allele copies in population
        count = RandomDist::binomial(n, p);
		
		// reduce residual offspring population size to fill
		n -= count;
	}
	
	assert(n == 0);
	assert(total == 0);
}


std::vector<unsigned int> RandomDist::multinomialByValue(const std::vector<unsigned int>& pop, int n) {	
	std::vector<unsigned int> res = pop;
	
	RandomDist::multinomial(res, n);
	
	return res;
}