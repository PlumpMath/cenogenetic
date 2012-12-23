#include "chromosome.h"

namespace ceno {
using f_t  = ternary_func<double>;
using ntides_t = Nucleotides<f_t >;
using allele_t = Allele<f_t>;

ntides_t g;
template<> ntides_t* allele_t::ntides = &g;

}

using namespace ceno;

int main() {
	//if (argc < 2) return 0;
	g.term_names = {"x", "y", "z"};
	g.terms = { 2., 4.5, 7.6, .34 };
	g.funcs.push_back(f_t([](double a, double b, double c) { return a * b * c;}));
	g.func_names.push_back("mult");
	g.funcs.push_back(f_t([](double a, double b, double c) { return a + b + c;}));
	g.func_names.push_back("add");
	g.funcs.push_back(f_t([](double a, double b, double c) { return a - b - c;}));
	g.func_names.push_back("sub");
g.funcs.push_back(f_t([](double a, double b, double c) { return a > b ? a > c ? a : c : b > c ? b : c;}));
	g.func_names.push_back("max");
g.funcs.push_back(f_t([](double a, double b, double c) { return a < b ? a < c ? a : c : b < c ? b : c;}));
	g.func_names.push_back("min");
g.funcs.push_back(f_t([](double a, double b, double c) { return a > b ? a < c ? a : b < c ? b : c : a > c ? c : a;}));
	g.func_names.push_back("middle");
	Chromosome<f_t> als;
	for (int i = 1; i < 5; i++) {
		als.clear();
		allele_t::build(back_inserter(als), CREATE_FULL, i);
		std::cout << "FULL:\n\n";
		std::cout << sexp << als << std::endl << std::endl;
		//std::cout << postfix << als << std::endl << std::endl;
		als.clear();
		allele_t::build(back_inserter(als), CREATE_GROW, i);
		std::cout << "\nGROW:\n\n";
		std::cout << sexp << als << std::endl << std::endl;
		//std::cout << postfix << als << std::endl << std::endl;
	}
}
