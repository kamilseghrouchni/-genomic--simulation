#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <vector>


/*!
  This is a random number class based on standard c++-11 generators
 */
class RandomDist {
    
public:

	/*!
	  Initializes the generator \ref rng with the Mersenne twister *mt19937* engine. 
	  Must provide mean *m*, standard deviation *sd* and sample size *ns*, 
	  set whether to use normal distribution (default uniform) and can provide a seed *s*.
	 */
    RandomDist(double m, double sd, int ns, bool n = false);
    

	/*!
	  Returns a vector of random doubles corresponding to the parameters set in the constructor.
	*/
    std::vector<double> generate_numbers();


	/** \brief Get a number following a binomial distribution
	 *
	 * Uses a mt19937 mersenne twister.
	 *
	 * */
    static int binomial(int n, double p);

    
    /** \brief Get an integer following a uniform distribution in the range [min, max]
	 *
	 * \param min		minimal value of distribution
	 * \param max		maximal value of distribution
	 * 
	 * \return Random integer in the given range
	 * */
    static int uniformIntSingle(int min, int max);


    /** \brief Get a double following a uniform distribution in the range [min, max]
	 *
	 * \param min		minimal value of distribution
	 * \param max		maximal value of distribution
	 * 
	 * \return Random double in the given range
	 * */
    static double uniformDoubleSingle(double min, double max);
    

    /** \brief Fill a vector with integers following a uniform distribution in the range [min, max]
	 *
	 * \param toFill	vector to fill, presupposedly initialized to a certain size
	 * \param min		minimal value of distribution
	 * \param max		maximal value of distribution
	 * */
    static void uniformIntVector(std::vector<int>& toFill, int min, int max);


    /** \brief Fill a vector with doubles following a uniform distribution in the range [min, max]
	 *
	 * \param toFill	vector to fill, presupposedly initialized to a certain size
	 * \param min		minimal value of distribution
	 * \param max		maximal value of distribution
	 * */
    static void uniformDoubleVector(std::vector<double>& toFill, double min, double max);

    
    /** \brief Generate a offspring population based on the parent population
     * 
     * Replaces parent population
     * Calls multinomial(pop, {size of population});
	 *
	 * \param pop		parent population
	 * */
    static void multinomial(std::vector<unsigned int>& pop);


    /** \brief Generate a offspring population based on the parent population
     * 
     * Replaces parent population
     * Calls multinomial(pop, {size of population});
	 *
	 * \param pop		parent population
	 * \param n			the size of the child population
	 * */
   

    static void multinomial(std::vector<unsigned int>& pop, int n);
    /** \brief Generate a offspring population based on the parent population
     * 
     * Generates new vector, does not overwrite passed argument
     * Calls multinomial(pop, {size of population});
	 *
	 * \param pop		parent population
	 * \param n			the size of the child population
	 * */
    static std::vector<unsigned int> multinomialByValue(const std::vector<unsigned int>& pop, int n);
     
private:

	//!< Random device
	static std::random_device rd;


	//!< Random number generator
	static std::mt19937 rng;


	//!< Fill a vector with uniform values
    void uniform(std::vector< double >&);


    //!< Fill a vector with normally distributed values
    void normal(std::vector< double >&);
    

    //!< Mean of the distribution
    double mean;


    //!< Standard deviation of the distribution
    double sd;


    //!< Size of the sample
    int nsample;


    //!< Flag to indicate normal oder uniform distribution (true for normal)
    bool normdist;
};

#endif
