#ifndef SMITH_WATERMAN_HPP
#define SMITH_WATERMAN_HPP

namespace alignment {

	double smithwaterman(std::string&, std::string&, double match = 2.0, double mismatch = -0.5, double indel = -1.0, bool verbose = false);

}

#endif
