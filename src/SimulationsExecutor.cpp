#include <algorithm>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <thread>
#include <ctime>
#include "SimulationsExecutor.hpp"
#include "Random.hpp"


SimulationsExecutor::SimulationsExecutor(std::string input, std::string fasta)
  : data(input, fasta)
{
	int allelesCountSum = 0;
	for (auto& alleleCount : data.getAllelesCount())
		allelesCountSum += alleleCount;

	assert(data.getPopulationSize() == allelesCountSum);

	switch (data.getExecutionMode()) {
		case _EXECUTION_MODE_MUTATIONS_:
			// generate nucleotide mutation probabilities according to model
			generateMutationRates();
			break;

		case _EXECUTION_MODE_MIGRATION_:
			generateSubPopulations();
			break;
			
		default:
			break;
	}
	    
    // init number of threads
    nThreads = std::thread::hardware_concurrency();
    if (nThreads == 0) {
		nThreads = 4;
		std::cout << "No hardware info detected, running on 4 threads" << std::endl;
	} else {
		std::cout << "Running on " << nThreads << " threads by system recommandation" << std::endl;
	}
	
	// output vals
	outputVals = std::vector< std::vector<std::string> >(
            (unsigned long) (data.getNbGenerations() + 2),
			std::vector<std::string>((unsigned long) data.getNbReplicates())
	);
}


void SimulationsExecutor::execute() {
	// chrono
	time_t t1 = time(0);
	
	// create simulation partition
	std::vector<int> nbSimulations(nThreads, data.getNbReplicates() / nThreads);
	
	// assign remaining simulations
	int rest = data.getNbReplicates() % nThreads;
	while (rest > 0) {
		++nbSimulations[rest];
		--rest;
	}
	
	// create correct number of threads
	std::vector<std::thread> threads(nThreads);

	// init each thread
	int minSimulationIdx = 0;
	for (size_t i = 0; i < nThreads; ++i) {
		threads[i] = std::thread([=] {
			runSimulation(nbSimulations[i], minSimulationIdx);
		});
		
		minSimulationIdx += nbSimulations[i];
	}

	// join threads (wait for every thread to end before ending main thread)
	for (auto& th : threads) th.join();
	
	// write data to ouput file
	writeData();

	// get end of the simulation
	time_t t2 = time(0);

	std::cout << "[done in " << t2 - t1 << " s]" << std::endl;
}


Simulation SimulationsExecutor::createSimulation() const {
	switch (data.getExecutionMode()) {
		case _EXECUTION_MODE_MUTATIONS_:
			return Simulation(data.getAlleles(), data.getAllelesCount(), data.getMutationRates(), nuclMutationProbs);

		case _EXECUTION_MODE_MIGRATION_:
			return Simulation(data.getAlleles(), subPopulations, migrationRates, data.getIsDetailedOutput());
			
		case _EXECUTION_MODE_SELECTION_:
			return Simulation(data.getAlleles(), data.getAllelesCount(), data.getSelections());

		case _EXECUTION_MODE_BOTTLENECK_:
			return Simulation(data.getAlleles(), data.getAllelesCount(), 
				data.getBottleneckStart(), data.getBottleneckEnd(), data.getPopReduction());

		case _EXECUTION_MODE_NONE_:
		default:
			return Simulation(data.getAlleles(), data.getAllelesCount());
	}
}


void SimulationsExecutor::runSimulation(int nSimulations, int firstSimulationIdx) {
	if (nSimulations == 0) return;
	
	// generate container for states of simulation
	int T = data.getNbGenerations();
	
	for (int i = firstSimulationIdx; i < nSimulations + firstSimulationIdx; ++i) {
		// create new simulation
		Simulation simul = createSimulation();		

		// write initial allele frequencies
		outputVals[0][i] = simul.getAlleleFqsForOutput();

		int t = 0;
		while (t < T) {
			// update simulation
			simul.update(t);

			// increment clock
			++t;

			// write allele frequencies
			outputVals[t][i] = simul.getAlleleFqsForOutput();
		}

		// write final line: allele identifiers
		outputVals[t + 1][i] = simul.getAlleleStrings();

		// properly format output (add blanks if there were new mutations
		if (data.getExecutionMode() == _EXECUTION_MODE_MUTATIONS_) {
			size_t lineLength = simul.getAlleleStrings().size();
			size_t precision = simul.getPrecision();

			if (outputVals.front()[i].size() != outputVals.back()[i].size()) { // check for any mutations
				for (size_t j = 0; j < outputVals.size(); ++j) {
					auto& state = outputVals[j][i];
					
					if (state.size() < lineLength) {
						std::stringstream ss;
						ss << state;

						// add additional zeroes to end of line
						while (ss.tellp() < (int) lineLength) {
							ss << _OUTPUT_SEPARATOR_ << std::setprecision((int) precision) << std::fixed << 0.0;
						}

						state = ss.str();
					}
				}
			}
		}
	}
}


void SimulationsExecutor::writeData() {
	// open result file
    results.open("results.txt");
    
	// write to result file
	for (int i = 0; i < (int) outputVals.size(); ++i) {
		writeAlleleFqs(i, outputVals[i]);	
	}
}


void SimulationsExecutor::writeAlleleFqs(int step, const std::vector<std::string>& alleleFqs) {
	results << step;
	if (data.getNbGenerations() > 998 && step < 1000) {
		if (step < 10)
			results << " ";
		if (step < 100)
			results << " ";
		
		results << " ";
	}
		
	results << '\t';
	
	for (auto const& data : alleleFqs) {
		results << data << '\t';
	}

	results << '\n';
}


void SimulationsExecutor::generateMutationRates() {
	switch (data.getMutationModel()) {
		case _MUTATION_MODEL_CANTOR_:
			{
				//std::cout << "Cantor model" << std::endl;

				double p = 1.0 / 3.0;

				nuclMutationProbs = { {
					{ { 0.0, p, p, p } },
					{ { p, 0.0, p, p } },
					{ { p, p, 0.0, p } },
					{ { p, p, p, 0.0 } }
				} };
			}
			break;

		case _MUTATION_MODEL_KIMURA_:
			{
				//std::cout << "Kimura model" << std::endl;

				double transition = data.getKimuraDelta();
				double transversion = (1.0 - transition) / 2.0;

				nuclMutationProbs = { {
					//  A    C			   G		   T
					{ { 0.0, transversion, transition, transversion } },
					{ { transversion, 0.0, transversion, transition } },
					{ { transition, transversion, 0.0, transversion } },
					{ { transversion, transition, transversion, 0.0 } }
				} };
			}
			break;

		case _MUTATION_MODEL_FELSENSTEIN_:
			{
				//std::cout << "Felsenstein model" << std::endl;

				std::vector<double> consts = data.getFelsensteinConstants();

				for (auto& c : consts) {
					assert(c != 1.0);
				}

				double pA = consts[Nucl::Nucleotide::A] / (1.0 - consts[Nucl::Nucleotide::A]);
				double pC = consts[Nucl::Nucleotide::C] / (1.0 - consts[Nucl::Nucleotide::C]);
				double pG = consts[Nucl::Nucleotide::G] / (1.0 - consts[Nucl::Nucleotide::G]);
				double pT = consts[Nucl::Nucleotide::T] / (1.0 - consts[Nucl::Nucleotide::T]);

				nuclMutationProbs = { {
					{ { 0.0, pC, pG, pT } },
					{ { pA, 0.0, pG, pT } },
					{ { pA, pC, 0.0, pT } },
					{ { pA, pC, pG, 0.0 } }
				} };

				// normalize table rows
				for (std::size_t i = 0; i < nuclMutationProbs.size(); ++i) {
					double sum = 0.0;
					for (std::size_t j = 0; j < nuclMutationProbs[i].size(); ++j) {
						sum += nuclMutationProbs[i][j];
					}

					assert(sum > 0.0);

					for (std::size_t j = 0; j < nuclMutationProbs[i].size(); ++j) {
						nuclMutationProbs[i][j] /= sum;
					}
				}
			}
			break;

		default:
			//std::cout << "No mutation model" << std::endl;
			break;
	}

	// display the prob matrix
	for (auto& pX : nuclMutationProbs) {
		for (auto& p : pX) {
			std::cout << p << '\t';
		}
		std::cout << std::endl;
	}
}


void SimulationsExecutor::generateSubPopulations() {
	// subpopulations generation
	subPopulations = std::vector< std::vector<unsigned int> >(data.getNbAlleles(), 
															  std::vector<unsigned int>(data.getNbAlleles(), 0));
		
    for (auto it = data.getAllelesCount().begin(); it != data.getAllelesCount().end(); ++it) {
        long int idx = it - data.getAllelesCount().begin();

        subPopulations[idx][idx] += rint(*it);

        assert(subPopulations[idx][idx] > 0);
    }

    // display the subpopulations
   /* std::cout << "Initial subpopulations" << std::endl;
	for (auto& pop : subPopulations) {
		for (auto& count : pop) {
			std::cout << count << '\t';
		}
		std::cout << std::endl;
	}
	std::cout << "---------" << std::endl;*/


    // migration rate generation
    migrationRates = std::vector<std::vector<unsigned int> >(data.getNbAlleles(),
					 std::vector<unsigned int>(data.getNbAlleles(), 0));

    starCenter = std::rand() % data.getNbAlleles();
	assert(starCenter <= data.getNbAlleles());
	assert(starCenter >= 0);

    for (size_t i(0); i < migrationRates.size() - 1; ++i) {

		size_t rate = 0;

		switch (data.getMigrationMode()) {
			case _MIGRATION_MODE_INPUT_USER_:
				{
					// check if the migration rate index is valid
					if (i < data.getMigrationRates().size()) {
						rate = (size_t) data.getMigrationRates()[i];
					}
				}
				break;

			case _MIGRATION_MODE_RANDOM_:
				{
					// take the value of the smallest subgroup still to distribute
			    	size_t maxMoving(data.getAllelesCount().back());
					
					std::vector<unsigned int> remaining(data.getAllelesCount().begin() + i, data.getAllelesCount().end());
					for (auto& elt : remaining) {
						maxMoving = (size_t ) std::min(std::max((int) elt, 0), (int) maxMoving);
					}

					maxMoving = (size_t) (maxMoving * 1.0 / (migrationRates.size() - 1 - i)); // will never be a division by 0 

					// randomly chosen rate
		            if (maxMoving > 1)
		                rate = (size_t) RandomDist::uniformIntSingle(1,(int) maxMoving - 1);
		            else
		                rate = maxMoving;
	            }
	            break;

	        default:
	         	break;
		}

        for (size_t j = i + 1; j < migrationRates[i].size(); ++j) {
			switch (data.getMigrationModel()) {
				case _MIGRATION_MODEL_COMPLETE_GRAPH_:
					{
						//exchanges between all and every subpopulation
						migrationRates[i][j] = (unsigned int) rate;
						migrationRates[j][i] = (unsigned int) rate;
					}
					break;

				case _MIGRATION_MODEL_STAR_:
					{
						// only the subpopulation in the center exchanges with all the others
						if (starCenter == i) {
							migrationRates[i][j] = (unsigned int) rate;
							migrationRates[j][i] = (unsigned int) rate;
						}
					}
					break;

				case _MIGRATION_MODEL_RING_:
					{
						// exchanges between "neighbor subpopulations"
						if (j == i + 1 || (j == migrationRates[i].size() - 1 && i == 0)) {
							migrationRates[i][j] = (unsigned int) rate;
							migrationRates[j][i] = (unsigned int) rate;
						}
					}
				break;

				default:
					break;
			}
		}
	}

	// check the migration rate table
	for (size_t i = 0; i < subPopulations.size(); ++i) {
		int popSum = 0;
		for (size_t j = 0; j < subPopulations[i].size(); ++j) {
			popSum += subPopulations[i][j];
		}

		int migSum = 0;
		for (size_t j = 0; j < migrationRates[i].size(); ++j) {
			migSum += migrationRates[i][j];
		}

		if (migSum > popSum) {
			// there are too many outgoing individuals!
			// we need to reduce it to the maximal possible value

			size_t j = 0;
			while (migSum > popSum && migSum > 0) {
				if (migrationRates[i][j] > 0) {
					--migrationRates[i][j];
					--migSum;
				}

				++j;
				j %= migrationRates[i].size();
			}
		}
	}

	// display the mutation rates
	/*std::cout << "Migration rate table" << std::endl;
	for (auto& mig : migrationRates) {
		for (auto& m : mig) {
			std::cout << m << '\t';
		}
		std::cout << std::endl;
	}*/
}


size_t SimulationsExecutor::getStarCenter() const {
	return starCenter;
}
