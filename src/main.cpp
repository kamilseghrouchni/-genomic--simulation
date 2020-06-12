#include "SimulationsExecutor.hpp"
#include "Data.hpp"


int main(int argc, char** argv) {
	std::string inputFileName = argc > 1 ? argv[1] : "../data/input.txt";	
	std::string fastaFileName = argc > 2 ? argv[2] : "";
	
	SimulationsExecutor simulationsExecutor(inputFileName, fastaFileName);
	simulationsExecutor.execute();
	
	return 0;
}
