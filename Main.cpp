#include<iostream>
#include "Header.h"
//#include "Measure.h"
using namespace Simplification;
//using namespace Variable;
extern std::string inputcloud;
int main() {
	//inputcloud = "S0001A0022_XYZ";
	std::cout << "input file:";
	std::cin >> inputcloud;
	int func = 0;
	std::cout << "choose function:" << std::endl;
	std::cout << "1.Visualize\n2.calculation the curvature\n3.simplification" << std::endl;
	std::cin >> func;
	switch (func)
	{
	case 1:
		visualization Vs;
		Vs.Visualize();
		break;

	case 2:
		Curvature Cr;
		Cr.mean_curv();
		break;

	case 3:
		simplicate Sp;
		float alpha;
		std::string Out;
		std::cout << "enter your alpha(調整簡化程度):" << std::endl;
		std::cin >> alpha;
		std::cout << "enter your output filename:" << std::endl;
		std::cin >> Out;
		Sp.simp(alpha,Out);
		break;
	}
	
	
	

	
	
	return 0;
}