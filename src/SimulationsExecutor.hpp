#ifndef SIMULATIONS_EXECUTOR_H
#define SIMULATIONS_EXECUTOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <mutex>
#include "Simulation.hpp"
#include "Data.hpp"
#include "Globals.hpp"


/** \brief Class representing a SimulationsExecutor
 * 
 * This class is a wrapper for the execution of multiple single
 * simulations with the same initial parameters (e.g. to generate statistics).
 * The Simulations are run on different threads and their outputs stored
 * in a common file.  
 * */
class SimulationsExecutor {
	
public:
	
	/** \brief SimulationsExecutor constructor
	 * 
	 * Initialises a new SimulationsExecutor, a wrapper for a series of Simulations
	 * uses the given parameters to construct a Data object to read the simulation params.
	 *
	 * \param input 		the path of the input file to be read
	 * \param fasta 		the path of the fasta file to be read
	 * */
	SimulationsExecutor(std::string input, std::string fasta);
	

	//!< Because of the threads, we do not allow any copies
	SimulationsExecutor(const SimulationsExecutor& other) = delete;
	

	//!< Because of the threads, we do not allow any copies
	SimulationsExecutor& operator=(const SimulationsExecutor& other) = delete;
	

	/** \brief Start the execution of the simulations
	 * 
	 * Invoking this method will create the threads and execute them.
	 * */
	void execute();

protected:

	/** \brief Generate a new Simulation based on the given parameters
	 *
	 * \return A new Simulation based on the user's paramters
	 * */
	Simulation createSimulation() const;
	

	/** \brief Run a simulation
	 * 
	 * This method is executed by a thread.
	 * 
	 * \param nSimulations			number of simulations to be run
	 * \param firstSimulationIdx	simulation index offset (relevant for output)
	 * */
	void runSimulation(int nSimulations, int firstSimulationIdx);
	

	/** \brief Write data to the result file
	 * */
	void writeData();
	

	/** \brief Write one step of all simulations to the result file
	 * 
	 * Wrties the data passed as argument in the result file.
	 * 
	 * \param step			the step number of the simulation to be written
	 * \param alleleFqs		a vector of strings (each being a formatted list of allele frequencies)
	 * */
	void writeAlleleFqs(int step, const std::vector<std::string>& alleleFqs);
	

	/** \brief Generate table for nucleotide mutation rates based on user input data
	 * 
	 * Automatically selects the correct model based on user data
	 * */
	void generateMutationRates();


	 /** \brief Generate the subpopulations for a Simulation with migration
    * 
	 * Automatically selects the correct model based on user data
	 * */
	void generateSubPopulations();


	 /** \brief Get the migration model used, if any
	 *
	 * \return An integer representing the migration model
	 * */
	int getMigrationModel() const ;


	/** \brief Get the index of the population at the center of a star-shaped migration pattern
	 *
	 * \return An size_t representing the index of the population at the center of a star-shaped migration pattern
	 * */
	size_t getStarCenter() const;	

private:

	//!< Data object containing all user params
	Data data;
	

	//!< Table of mutation probabilities
	std::array< std::array<double, Nucl::Nucleotide::N >, Nucl::Nucleotide::N > nuclMutationProbs;
	
	
	//!< Table containing the sub-populations
    std::vector< std::vector<unsigned int> > subPopulations;

    
    //!< Table containing migration rates for each sub group
    std::vector< std::vector<unsigned int> > migrationRates;


    //!< The index of the population at the center of a star-shaped migration pattern
    size_t starCenter;
	

	//!< Result file
	std::ofstream results;
	
	//!< Number of threads
	unsigned int nThreads;
	
	//!< Output values
	std::vector< std::vector<std::string> > outputVals;
};

#endif
