#pragma once
//Header.h
#ifndef HEADER_H
#define HEADER_H

#include <string>
#include<vector>

namespace Simplification {

	static class visualization {
	public:
		std::string Visualize();
	};

	static class Curvature {
	public:
		float mean_curv();
	};

	static class simplicate {
	public:
		std::string simp(float alpha, std::string Out);
	};
}

#endif HEADER_H