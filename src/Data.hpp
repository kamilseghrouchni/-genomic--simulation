#ifndef DATA_H
#define DATA_H

#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <functional> 
#include "Globals.hpp"


/** \brief Class regrouping the data necessary to run a simulation
 *
 * Collecting the data requires : a user input file (.txt) and a fasta file (.fa)
 *
 * */
class Data {

public:

	/** \brief Data constructor
	 *
	 * Initialises the name of the data files, according to the parameters
	 * Loads all the parameters form the file
	 *
	 * \param input 		the path of the input file to be read
	 * \param fasta 		the path of the fasta file to be read
	 * */
	Data(std::string input, std::string fasta);
	
	
	/** \brief Getter of the size of the populationSize
	 *
	 * 	\return populationSize, an int
	 * */
	int getPopulationSize() const;


	/** \brief Getter of the number of generations
	 *
	 * 	\return numberGenerations, an int
	 * */
	int getNbGenerations() const;


	/** \brief Getter of the number of replicates
	 *
	 * 	\return replicates, an int
	 * */
	int getNbReplicates() const;


	/** \brief Getter of the number of alleles
	 *
	 * 	\return numberAlleles, a size_t
	 * */
	size_t getNbAlleles() const;


	/** \brief Getter of the vector of the allele counts
	 * 
	 * 	\return a vector of unsigned integers, the allele counts
	 * */
	const std::vector<unsigned int>& getAllelesCount() const;


	/** \brief Getter of the vector of the alleles
	 * 
	 * 	\return a vector of strings, the alleles
	 * */
	const std::vector<std::string>& getAlleles() const;


	/** \brief Get the marker sites on the alleles
	 * 
	 * \return a vector of unsigned integers, the marker sites
	 * */
	const std::vector<unsigned int>& getMarkerSites() const;

	
	/** \brief Get the execution mode of a Simulation (mutation, migration, ...)
	 * 
	 * \return An int whose meaning is defined in Globals.hpp
	 * */
	int getExecutionMode() const;


	/** \brief Getter of the sites mutations probabilities
	 *
	 * 	\return mutations, a vector of double
	 * */
	const std::vector<double>& getMutationRates() const;


	/** \brief Get the mutation model to use for a Simulation (Cantor, Kimura, ...)
	 * 
	 * */
	int getMutationModel() const;

	
	/** \brief Get the value necessary to create a Kimura mutation model
	 * 
	 * \return The value of delta for a Kimura mutation model
	 * */
	double getKimuraDelta() const;
	

	/** \brief Get the values necessary to create a Felsenstein mutation model
	 * 
	 * \return The values of the constants for a Felsenstein model
	 * */
	const std::vector<double>& getFelsensteinConstants() const;


	/** \brief Get the migration model to use for a Simulation (comppleteGraph, Star, ...)
	 *
	 * */
    int getMigrationModel() const;


	/** \brief Get the migration mode to use for a Simulation (random, user input , ...)
	 *
	 * */
	int getMigrationMode() const;
	

	/** \brief Getter of the sites migration rates for each allele
     *
     * The user is supposed to know the number of different allele.
     * Each rate correspond to the amount of outgoing individuals.
     *
     * \return The migration rates, a vector of integers
     * */
    const std::vector<int>& getMigrationRates() const;
    

    /** \brief Get whether the output should be detailed if in migration mode
     * 
     * */
    bool getIsDetailedOutput() const;


	/** \brief Get selection rates
	 * */
	const std::vector<double>& getSelections() const;
	

	/** \brief Getter of the reduction of the population size factor during the bottleneck
	 *
	 * 	\return popReduction, a double
	 * */
	double getPopReduction() const;

	 
	/** \brief Getter of the start of the bottleneck
	 *
	 * 	\return bottleneckStart, an int
	 * */
	int getBottleneckStart() const;
	 

	/** \brief Getter of the end of the bottleneck
	 *
	 * 	\return bottleneckEnd, an int
	 * */
	int getBottleneckEnd() const;
	
	
	/** \brief Utility function to transfrom strings to integers for use in switch statements
	 * 
	 * */
	static constexpr unsigned int str2int(const char* str, int h = 0) {
		return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
	} 

    
    /** \brief Extracts a value from a line
	 * 
	 * \param into			the variable to store the extracted value into
	 * \param line			the line to read from, a string
	 * \param extractFn		the function to call to cast the value contained in the string to the corresponding arithmetic type
	 * */
    template<typename T>
	void extractValue(T& into, std::string line, std::function< T (const std::string&) > extractFn) {
		std::string key, strValue;

		std::stringstream ss(line);
		std::getline(ss, key, _INPUT_DECLARATION_); // separate key from value
		std::getline(ss, strValue); // store value in strValue
		
		try {
			// cast value
			into = extractFn(strValue);
			
		} catch (std::invalid_argument& e) {	
			std::cerr << e.what() << std::endl;
			throw e.what();
		}
	}
	

	/** \brief Extracts a vector of values from a line
	 * 
	 * \param into			the variable to store the extracted values into
	 * \param line			the line to read from, a string
	 * \param extractFn		the function to call to cast the value contained in the string to the corresponding arithmetic type
	 * */
	template<typename T>
	void extractValues(std::vector<T>& into, std::string line, std::function< T (const std::string&) > extractFn) {
		std::string key, strValue;

		std::stringstream ss(line);
		std::getline(ss, key, _INPUT_DECLARATION_); // separate key from value

		while (std::getline(ss, strValue, _INPUT_SEPARATOR_)) {
			try {
				//add elements found between each separators
				T extractedValue = extractFn(strValue);
				into.push_back(extractedValue);

			} catch (std::invalid_argument& e) {
				std::cerr << e.what() << std::endl;
				throw e.what();
			}
		}
	}

protected:

	/** \brief Collects all the data required for the program
	 *
	 * Collects inputs from the user and the fasta files
	 * */
	void collectAll();


	/** \brief Collects data from the user file
	 *
	 * Reads the number of generations, the marker sites, the number of getReplicates
	 * Also reads the migration probabilities of each nucleotide
	 * */
	void collectUserFile(std::ifstream& file);


	/** \brief Checks the data from the user file
	 *
	 * Reviews the read values and evalutes if the parameters are correct to run a Simulation
	 * */
	void checkUserFile();


	/** \brief Collects data from the fasta file
	 *
	 * Calculates the number of individuals/size of the population
	 * Registers the allele sequences of all the individuals
	 * Counts the number of different alleles and their frequencies
	 * Initialises the name of the input files
	 *
	 * */
	void collectFastaFile(std::ifstream& file);


	/** \brief Checks the data from the fasta file
	 *
	 * Reviews the read data
	 * Generates the alleles based on the read data
	 *
	 * */
	void checkFastaFile();

	
private:

	//!< Name of the user input file, a string
	std::string inputName;

	
	//!< Name of the fasta file, a string
	std::string fastaName;

	
	//!< Flag for operating mode
	bool withFasta;


	//!< Size of the population, an int
	int populationSize;

	
	//!< Number of generations (number of simulation steps), an int
	int nbGenerations;
	
	
	//!< Number of replicates of the simulation, an int
	int nbReplicates;

	
	//!< Vector containing the allele frequencies when running without fasta file
	std::vector<double> allelesFqs;

	
	//!< List of alleles in the Simulations
	std::vector<std::string> alleles;
	
	
	//!< Vector of double containing the allele frequencies of the fasta file
	std::vector<unsigned int> allelesCount;

	
	//!< Vector of double containing the user marker sites
	std::vector<unsigned int> markerSites;

	
	//!< List of strings containing the allele sequences of all the individuals of the simulation
	std::list<std::string> sequences;


	//!< Execution mode (param to use)
	int executionMode;
	

	//!< Vector of double containing the mutations probabilities of the marker sites
	std::vector<double> mutationRates;
	
	
	//!< Mutation model (simple, kimura, felsenstein)
	int mutationModel;
	
	
	//!< Kimura model
	double kimuraDelta;
	
	
	//!< Felsenstein model
	std::vector<double> felsensteinConstants;


	//!< Migration  model (complete graph, ring, star)
    int migrationModel;
    
    
    //!< Migration  mode (user input, random)
	int migrationMode;
    
    
    //!< Migration rates
    std::vector<int> migrationRates;
    
    
    //!< Flag for the detailed output or not
    bool isMigrationDetailedOutput;


	//!< Vector of double containing the selection probabilities of the alleles
	std::vector<double> selections;

	
	//!< Bottleneck population reduction factor
	double popReduction;
	
	
	//!< Bottleneck start time
	int bottleneckStart;
	
	
	//!< Bottleneck stop time
	int bottleneckEnd;
};

#endif
