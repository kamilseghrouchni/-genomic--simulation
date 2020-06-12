#ifndef GLOBALS_H
#define GLOBALS_H

#include <map>
namespace Nucl {
	typedef enum Nucleotide { A, C, G, T, N } Nucleotide;
	
	const std::map<char, Nucleotide> fromChar = { 
		{'a', A}, {'A', A},
		{'c', C}, {'C', C},
		{'g', G}, {'G', G},
		{'t', T}, {'T', T},
		{'n', N}, {'N', N}, {'x', N}, {'*', N}
	};
	
	const char toChar[6] = "ACGTN";
	const char possibleChars[5] = "ACGT";
}

#define _INPUT_COMMENT_ '#'
#define _INPUT_DECLARATION_ '='
#define _INPUT_SEPARATOR_ '|'

#define _FASTA_COMMENT_ '>'

#define _OUTPUT_SEPARATOR_ '|'
#define _MIN_OUTPUT_PRECISION_ 2

#define _ERROR_INPUT_UNREADABLE_CODE_ 1
#define _ERROR_INPUT_UNREADABLE_MSG_ "Input file impossible to open."
#define _ERROR_FASTA_UNREADABLE_CODE_ 2
#define _ERROR_FASTA_UNREADABLE_MSG_ "Fasta file impossible to open."
#define _ERROR_MUTATION_TARGET_UNFINDABLE_CODE_ 3
#define _ERROR_MUTATION_TARGET_UNFINDABLE_MSG_ "Did not find mutation target."
#define _ERROR_OUTPUT_BUFFER_CODE_ 4
#define _ERROR_OUTPUT_BUFFER_MSG_ "Error: trying to add data of a step that has already been written to the result file."

#define _ERROR_MARKER_SITE_OUT_OF_BOUNDS_CODE_ 5
#define _ERROR_MARKER_SITE_OUT_OF_BOUNDS_MSG_ "Error: at least one marker site is out of bounds. Please select values between 0 and the line length - 1."

#define _ERROR_POPULATION_SIZE_TOO_SMALL_CODE_ 6
#define _ERROR_POPULATION_SIZE_TOO_SMALL_MSG_ "Error: the population should have at least one individual."

#define _ERROR_NO_INITIAL_FREQUENCIES_CODE_ 7
#define _ERROR_NO_INITIAL_FREQUENCIES_MSG_ "Error: no initial frequencies were specified."

#define _ERROR_NO_EXECUTION_MODE_CODE_ 8
#define _ERROR_NO_EXECUTION_MODE_MSG_ "Error: invalid execution mode."

#define _ERROR_NO_GENERATIONS_CODE_ 9
#define _ERROR_NO_GENERATIONS_MSG_ "Error: number of generations must be > 0."
#define _ERROR_NO_REPLICATES_CODE_ 10
#define _ERROR_NO_REPLICATES_MSG_ "Error: number of replicates must be > 0."

#define _ERROR_TOO_MANY_FELSENSTEIN_CONSTS_CODE_ 11

#define _ERROR__CODE_ 
#define _ERROR__MSG_ ""

#define _ERROR__CODE_ 
#define _ERROR__MSG_ ""

#define _ERROR__CODE_ 
#define _ERROR__MSG_ ""

#define _ERROR__CODE_ 
#define _ERROR__MSG_ ""

#define _ERROR__CODE_ 
#define _ERROR__MSG_ ""

#define _INPUT_KEY_GENERATIONS_ "GEN"
#define _INPUT_KEY_REPLICAS_ "REP"
#define _INPUT_KEY_POPULATION_SIZE_ "POPSIZE"
#define _INPUT_KEY_INITIAL_FREQ_ "FREQ"
#define _INPUT_KEY_MARKER_SITES_ "SITES"
#define _INPUT_KEY_MODE_ "MODE"

#define _EXECUTION_MODE_NONE_ 0
#define _EXECUTION_MODE_MUTATIONS_ 1
#define _EXECUTION_MODE_MIGRATION_ 2
#define _EXECUTION_MODE_SELECTION_ 3
#define _EXECUTION_MODE_BOTTLENECK_ 4

#define _INPUT_KEY_MUTATION_RATES_ "MUT"
#define _INPUT_KEY_MUTATION_KIMURA_ "MUT_KIMURA"
#define _INPUT_KEY_MUTATION_FELSENSTEIN_ "MUT_FELSENSTEIN"

#define _INPUT_KEY_MIGRATION_MODEL_ "MIG_MODEL"
#define _INPUT_KEY_MIGRATION_RATES_ "MIG_RATES"
#define _INPUT_KEY_MIGRATION_DETAILED_OUTPUT_ "MIG_DETAILED_OUTPUT"

#define _INPUT_KEY_SELECTION_RATES_ "SEL"

#define _INPUT_KEY_BOTTLENECK_POPULATION_REDUCTION_ "POP_REDUCTION"
#define _INPUT_KEY_BOTTLENECK_START_TIME_ "POP_START"
#define _INPUT_KEY_BOTTLENECK_END_TIME_ "POP_END"


#define _DEFAULT_MUTATION_RATE_ 1E-6
#define _MUTATION_MODEL_NONE_ 0
#define _MUTATION_MODEL_CANTOR_ 1
#define _MUTATION_MODEL_KIMURA_ 2
#define _MUTATION_MODEL_FELSENSTEIN_ 3

#define _MIGRATION_MODEL_NONE_ 0
#define _MIGRATION_MODEL_COMPLETE_GRAPH_ 1
#define _MIGRATION_MODEL_STAR_ 2
#define _MIGRATION_MODEL_RING_ 3

#define _MIGRATION_MODE_NONE_ 0
#define _MIGRATION_MODE_INPUT_USER_ 1
#define _MIGRATION_MODE_RANDOM_ 2

#define _MIGRATION_OUTPUT_SEPARATOR_ "  "

#define strToInt [](const std::string& s) { return std::stoi(s); }
#define strToUnsignedInt [](const std::string& s) { return (unsigned int) std::stoi(s); }
#define strToDouble [](const std::string& s) { return std::stod(s); }

#endif
