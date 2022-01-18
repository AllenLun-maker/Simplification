#pragma once
//Header.h
#ifndef HEADER_H
#define HEADER_H

#include <string>

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
		std::string simp(float alpha);
	};
}

#endif HEADER_H