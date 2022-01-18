#include "Header.h"
//#include "Measure.h"
using namespace Simplification;
//using namespace Variable;
extern std::string inputcloud;
int main() {
	inputcloud = "S0001A0022_XYZ";
	visualization Vs;
	
	Curvature Cr;

	simplicate Sp;
	Sp.simp(0.005);
	return 0;
}