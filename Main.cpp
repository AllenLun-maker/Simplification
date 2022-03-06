#include<iostream>
#include "Def_Demo.cpp"
#include "Def_Cal.cpp"
//#include "Measure.h"
//using namespace Demo;

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
		Show demo;
		demo.Visualize();
		break;

	case 2:
		Calculation cal;
		cal.mean_curv();
		break;

	case 3:
		Calculation cal;
		float alpha;
		std::string Out;
		std::cout << "enter your alpha(調整簡化程度):" << std::endl;
		std::cin >> alpha;
		std::cout << "enter your output filename:" << std::endl;
		std::cin >> Out;
		cal.simp(alpha,Out);
		break;
	}
	
	return 0;
}