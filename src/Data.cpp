#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <algorithm>
#include "Data.hpp"
#include "Random.hpp"

using namespace std;


Data::Data(string input, string fasta)
  : inputName(input), fastaName(fasta), withFasta(fasta != ""),
	populationSize(0), nbGenerations(0),
	nbReplicates(0), executionMode(_EXECUTION_MODE_NONE_),
	mutationModel(_MUTATION_MODEL_NONE_), kimuraDelta(0.0),
	migrationModel(_MIGRATION_MODEL_NONE_), migrationMode(_MIGRATION_MODE_NONE_),
	isMigrationDetailedOutput(false),
	popReduction(1), bottleneckStart(0), bottleneckEnd(0)
{
	assert(allelesCount.empty());
	assert(markerSites.empty());
	assert(sequences.empty());
	assert(mutationRates.empty());
	assert(selections.empty());

	collectAll();
}


void Data::collectAll() {
	// read user file
	ifstream dataFile;
	dataFile.open(inputName, ifstream::in);

	if (dataFile.is_open()) {
		collectUserFile(dataFile);
	} else {
		cerr << _ERROR_INPUT_UNREADABLE_MSG_ << endl;
		exit(_ERROR_INPUT_UNREADABLE_CODE_);
	}

	// check the user file
	checkUserFile();


	if (withFasta) {
		// read fasta file
		ifstream fastaFile;
		fastaFile.open(fastaName, ifstream::in);

		if (fastaFile.is_open()) {
			collectFastaFile(fastaFile);
		} else {
			cerr << _ERROR_FASTA_UNREADABLE_MSG_ << endl;
			exit(_ERROR_FASTA_UNREADABLE_CODE_);
		}

		// check the fasta file
		checkFastaFile();
	}
}


void Data::collectUserFile(ifstream& file) {
	string line, key;

	while (getline(file, line)) {

		// remove whitespace
		line.erase(
			remove_if(line.begin(), line.end(), ::isspace),
			line.end()
		);

		// ignore comment lines
		if (line[0] == _INPUT_COMMENT_) continue;

		stringstream ss(line);
		getline(ss, key, _INPUT_DECLARATION_);

        switch (str2int(key.c_str())) {
			case str2int(_INPUT_KEY_GENERATIONS_):
                extractValue<int>(nbGenerations, line, strToInt);
                break;

			case str2int(_INPUT_KEY_REPLICAS_):
				extractValue<int>(nbReplicates, line, strToInt);
				break;

			case str2int(_INPUT_KEY_POPULATION_SIZE_):
				if (!withFasta)
					extractValue<int>(populationSize, line, strToInt);
				break;

			case str2int(_INPUT_KEY_INITIAL_FREQ_):
				if (!withFasta)
					extractValues<double>(allelesFqs, line, strToDouble);
				break;

            case str2int(_INPUT_KEY_MARKER_SITES_):
				extractValues<unsigned int>(markerSites, line, strToUnsignedInt);
				break;

			case str2int(_INPUT_KEY_MODE_):
				extractValue<int>(executionMode, line, strToInt);
				break;

			// MUTATIONS
			case str2int(_INPUT_KEY_MUTATION_RATES_):
				extractValues<double>(mutationRates, line, strToDouble);
                break;

            case str2int(_INPUT_KEY_MUTATION_KIMURA_):
				extractValue<double>(kimuraDelta, line, strToDouble);
				break;

			case str2int(_INPUT_KEY_MUTATION_FELSENSTEIN_):
				extractValues<double>(felsensteinConstants, line, strToDouble);
				break;

			// MIGRATIONS
			case str2int(_INPUT_KEY_MIGRATION_MODEL_):
				extractValue<int>(migrationModel, line, strToInt);
				break;

            case str2int(_INPUT_KEY_MIGRATION_RATES_):
                extractValues<int>(migrationRates, line, strToInt);
                break;

            case str2int(_INPUT_KEY_MIGRATION_DETAILED_OUTPUT_):
				{
					int detailedOutput = 0;
					extractValue<int>(detailedOutput, line, strToInt);

					isMigrationDetailedOutput = detailedOutput == 1;
				}
				break;

            // SELECTION
			case str2int(_INPUT_KEY_SELECTION_RATES_):
				extractValues<double>(selections, line, strToDouble);
				break;

			// BOTTLENECK
			case str2int(_INPUT_KEY_BOTTLENECK_POPULATION_REDUCTION_):
				extractValue<double>(popReduction, line, strToDouble);
				break;

			case str2int(_INPUT_KEY_BOTTLENECK_START_TIME_):
				extractValue<int>(bottleneckStart, line, strToInt);
				break;

			case str2int(_INPUT_KEY_BOTTLENECK_END_TIME_):
				extractValue<int>(bottleneckEnd, line, strToInt);
				break;

            default:
				break;
        }
	}
}


void Data::checkUserFile() {
	// check basic params
	if (!(nbGenerations > 0)) {
		cerr << _ERROR_NO_GENERATIONS_MSG_ << endl;
		cerr << "Running simulation with 10 generations" << endl;
		nbGenerations = 10;
	}

	if (!(nbReplicates > 0)) {
		cerr << _ERROR_NO_REPLICATES_MSG_ << endl;
		cerr << "Running 1 replicate" << endl;
		nbReplicates = 1;
	}

	if (!withFasta) {
		// check for sufficient population
		if (!(populationSize > 0)) {
			cerr << _ERROR_POPULATION_SIZE_TOO_SMALL_MSG_ << endl;
			exit(_ERROR_POPULATION_SIZE_TOO_SMALL_CODE_);
		}

		// check for alleles
		if (allelesFqs.empty()) {
			cerr << _ERROR_NO_INITIAL_FREQUENCIES_MSG_ << endl;
			exit(_ERROR_NO_INITIAL_FREQUENCIES_CODE_);
		}

		// generate allelesn and allele counts
		for (size_t i = 0; i < allelesFqs.size(); ++i) {
			string idx = to_string(i);

			alleles.push_back(idx);
			allelesCount.push_back((unsigned int) (allelesFqs[i] * populationSize));
		}

		// check for correct alleles count
		int allelesCountSum = 0;
		for (auto& alleleCount : allelesCount)
			allelesCountSum += alleleCount;

		if (allelesCountSum != populationSize) {
			cerr << "Error: number of individuals does not match population size. Is the sum of the allele frequencies 1?" << endl;
			exit(1);
		}
	}


	// set mutation model, if mutation_mode
	switch (executionMode) {
		case _EXECUTION_MODE_NONE_:
			break;

		case _EXECUTION_MODE_MUTATIONS_:
			if (!withFasta) {
				cerr << "Error: a simulation with mutations can only be done with a fasta file (alleles with real genotypes)." << endl;
				exit(0);
			}

			// default is cantor
			mutationModel = _MUTATION_MODEL_CANTOR_;

			if (kimuraDelta != 0.0) {
				if (kimuraDelta >= 1.0/3.0 && kimuraDelta <= 1.0) {
					mutationModel = _MUTATION_MODEL_KIMURA_;
				} else {
					cerr << "To use the Kimura model of mutations, the value of delta must be in the range [0.33333, 1]." << endl;
				}
			} else if (!felsensteinConstants.empty()) {
				size_t nConstants = felsensteinConstants.size();

				if (nConstants != (size_t) Nucl::Nucleotide::N) {
					cerr << "4 Felsenstein constants are required, but " << nConstants << " were provided." << endl;
					if (nConstants < (size_t) Nucl::Nucleotide::N) {
						cerr << "Filling in missing values uniformly." << endl;

						double sum = 0.0;
						for (auto& c : felsensteinConstants) sum += c;

						double fill = max((1.0 - sum) / ((int) Nucl::Nucleotide::N - nConstants), 0.0);
						for (size_t i = 0; i < (size_t) Nucl::Nucleotide::N - nConstants; ++i)
							felsensteinConstants.push_back(fill);

					} else {
						cerr << "Too many felsenstein constants: aborting." << endl;
						exit(_ERROR_TOO_MANY_FELSENSTEIN_CONSTS_CODE_);
					}
				}
				// check for correct constants
				double sum = 0.0;
				for (auto& c : felsensteinConstants) {

					// c can not be negative
					if (c < 0.0)
						cerr << "Felsenstein constants can not be negative. Taking the absolute value as constant." << endl;
					c = abs(c);

					// c can not be equal to 1
					sum += c;
				}

				// adjust terms
				if (sum < 1.0) {
					double adjust = (1.0 - sum) / felsensteinConstants.size();
					for (auto& c : felsensteinConstants)
						c += adjust;
				}

				// set mutation model
				if (!(sum > 1.0)) {
					mutationModel = _MUTATION_MODEL_FELSENSTEIN_;
				} else {
					cerr << "Error: sum of felstenstein constants is greater than 1.0." << endl;
				}
			}
			break;

		case _EXECUTION_MODE_MIGRATION_:
			migrationMode = migrationRates.empty() ? _MIGRATION_MODE_RANDOM_ : _MIGRATION_MODE_INPUT_USER_;
			break;

		case _EXECUTION_MODE_SELECTION_:
			break;

		case _EXECUTION_MODE_BOTTLENECK_:
			break;

		default:
			cerr << _ERROR_NO_EXECUTION_MODE_MSG_ << endl;
			exit( _ERROR_NO_EXECUTION_MODE_CODE_);
			break;
	}
}


void Data::collectFastaFile(ifstream& file) {
	populationSize = 0;

	string line;
	while (getline(file, line)) {

		if (line[0] == _FASTA_COMMENT_) {
			// increasing population thanks to the header
			++populationSize;
			continue;
		}

		string seq = "";
		size_t lineLength = line.size();
		for (auto& marker : markerSites) {
			// check that the marker is valid
			if (marker >= lineLength) {
				cerr << _ERROR_MARKER_SITE_OUT_OF_BOUNDS_MSG_ << endl;
				exit(_ERROR_MARKER_SITE_OUT_OF_BOUNDS_CODE_);
			}

			char c = line[marker];
			if (string(Nucl::possibleChars).find(c) != string::npos) {
				seq += c;
			}
			else {
				// if we have an unknown nucleotide, generate a valid one randomly
				seq += Nucl::toChar[RandomDist::uniformIntSingle(Nucl::Nucleotide::A, Nucl::Nucleotide::T)];
			}
		}

		sequences.push_back(seq);
	}
}


void Data::checkFastaFile() {
	// count the alleles
	list<string> uniqueSequences = sequences;

	// sort and remove the duplicates
	uniqueSequences.sort();
	uniqueSequences.unique();

	for (auto& seq : uniqueSequences) {
		auto count = (unsigned int) count_if(sequences.begin(), sequences.end(), [&](string allele) {
				return allele == seq;
		});

		allelesCount.push_back(count);
		alleles.push_back(seq);
	}
}


int Data::getPopulationSize() const {
	return populationSize;
}


int Data::getNbGenerations() const {
	return nbGenerations;
}


int Data::getNbReplicates() const {
	return nbReplicates;
}


size_t Data::getNbAlleles() const {
	return allelesCount.size();
}


const std::vector<unsigned int>& Data::getAllelesCount() const {
	return allelesCount;
}


const std::vector<std::string>& Data::getAlleles() const {
	return alleles;
}


const std::vector<unsigned int>& Data::getMarkerSites() const {
	return markerSites;
}


int Data::getExecutionMode() const {
	return executionMode;
}


const vector<double>& Data::getMutationRates() const {
	return mutationRates;
}


int Data::getMutationModel() const {
	return mutationModel;
}


double Data::getKimuraDelta() const {
	return kimuraDelta;
}


const std::vector<double>& Data::getFelsensteinConstants() const {
	return felsensteinConstants;
}


int Data::getMigrationModel() const {
	return migrationModel;
}


int Data::getMigrationMode() const {
	return migrationMode;
}


const std::vector<int>& Data::getMigrationRates() const {
    return migrationRates;
}


bool Data::getIsDetailedOutput() const {
	return isMigrationDetailedOutput;
}


const std::vector<double>& Data::getSelections() const {
	return selections;
}


double Data::getPopReduction() const {
	return popReduction;
}


int Data::getBottleneckStart() const {
	return bottleneckStart;
}


int Data::getBottleneckEnd() const {
	return bottleneckEnd;
}
