#include "chromosome.h"

namespace ceno {
using T = double;
using U = unsigned char;
using ufunc = unary_func<T>;
using bfunc = binary_func<T>;
using tfunc = ternary_func<T>;
//using f_t  = typename tfunc::f_t;
using ntides_t = Nucleotides<T>;
using allele_t = Allele<T, U>;
using chromosome_t = Chromosome<T, U>;
using genome_t = Genome<T, U>;
ntides_t g;
template<> ntides_t* allele_t::ntides = &g;

}

using namespace ceno;

int main(int argc, char**) {
	//no output w/o argument
	if (argc < 2) {
		std::cout.setstate(std::ios::failbit);
	}
	g.term_names = {"x", "y", "z"};
	g.terms = { 2., 4.5, 7.6, .34 };
	g.funcs.push_back(std::unique_ptr<bfunc>(new bfunc([](T a, T b) {
		return a * b;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](T a, T b, T c) {
		return a + b + c;
	})));
	g.funcs.push_back(std::unique_ptr<bfunc>(new bfunc([](T a, T b) {
		return a - b;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](T a, T b, T c) {
	return a > b ? a > c ? a : c : b > c ? b : c;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](T a, T b, T c) {
	return a < b ? a < c ? a : c : b < c ? b : c;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](T a, T b, T c) {
	return a > b ? a < c ? a : b < c ? b : c : a > c ? c : a;
	})));
	g.funcs.push_back(std::unique_ptr<ufunc>(new ufunc([](T a) {
		return a * a;
	})));
	g.func_names.push_back("mult");
	g.func_names.push_back("add");
	g.func_names.push_back("sub");
	g.func_names.push_back("max");
	g.func_names.push_back("min");
	g.func_names.push_back("middle");
	g.func_names.push_back("square");
	genome_t grow, full;
	for (int i = 1; i < 5; i++) {
		full.clear();
		allele_t::build(back_inserter(full), CREATE_FULL, i);
		std::cout << "Depth: " << i << "\nFULL:\n\n";
		std::cout << sexp << full << std::endl << std::endl;
		std::cout << postfix << full << std::endl << std::endl;
		std::cout << "evals to " << allele_t::eval(full.begin()) << std::endl;
		grow.clear();
		allele_t::build(back_inserter(grow), CREATE_GROW, i);
		std::cout << "\nGROW:\n\n";
		std::cout << sexp << grow << std::endl << std::endl;
		std::cout << postfix << grow << std::endl << std::endl;
		std::cout << "evals to " << allele_t::eval(grow.begin()) << std::endl;
		genome_t jr, sis;
		breed(full, grow, std::back_inserter(jr), std::back_inserter(sis), .5);
		std::cout << sexp << "OFFSPRING: " << jr << std::endl << sis << std::endl;
	}
}
