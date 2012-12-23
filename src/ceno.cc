#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>

#include "rand.h"
#include "chromosome.h"

namespace { /*anon*/
std::default_random_engine rand_engine(std::chrono::system_clock::to_time_t(
		std::chrono::system_clock::now()));
std::uniform_int_distribution<size_t> rand_distro;
}

namespace ceno {

size_t bounded_rand(size_t boundary) {
	return rand_distro(rand_engine, std::uniform_int_distribution<size_t>::param_type(0, boundary - 1));
}
std::ios_base& prefix(std::ios_base& str) {
	str.setf(std::ios_base::left, std::ios_base::adjustfield);
	return str;
}
std::ios_base& postfix(std::ios_base& str) {
	str.setf(std::ios_base::right, std::ios_base::adjustfield);
	return str;
}
std::ios_base& sexp(std::ios_base& str) { return prefix(str); }
std::ios_base& forth(std::ios_base& str) { return postfix(str); }
}
