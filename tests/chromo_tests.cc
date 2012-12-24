#include "chromosome.h"

namespace ceno {
using ufunc = unary_func<double>;
using bfunc = binary_func<double>;
using tfunc = ternary_func<double>;
using f_t  = typename tfunc::f_t;
using ntides_t = Nucleotides<double>;
using allele_t = Allele<double>;

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
	g.funcs.push_back(std::unique_ptr<bfunc>(new bfunc([](double a, double b) {
		return a * b;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](double a, double b, double c) {
		return a + b + c;
	})));
	g.funcs.push_back(std::unique_ptr<bfunc>(new bfunc([](double a, double b) {
		return a - b;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](double a, double b, double c) {
	return a > b ? a > c ? a : c : b > c ? b : c;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](double a, double b, double c) {
	return a < b ? a < c ? a : c : b < c ? b : c;
	})));
	g.funcs.push_back(std::unique_ptr<tfunc>(new tfunc([](double a, double b, double c) {
	return a > b ? a < c ? a : b < c ? b : c : a > c ? c : a;
	})));
	g.funcs.push_back(std::unique_ptr<ufunc>(new ufunc([](double a) {
		return a * a;
	})));
	g.func_names.push_back("mult");
	g.func_names.push_back("add");
	g.func_names.push_back("sub");
	g.func_names.push_back("max");
	g.func_names.push_back("min");
	g.func_names.push_back("middle");
	g.func_names.push_back("square");
	Chromosome<double> algrow,alfull;
	for (int i = 1; i < 5; i++) {
		alfull.clear();
		allele_t::build(back_inserter(alfull), CREATE_FULL, i);
		std::cout << "Depth: " << i << "\nFULL:\n\n";
		std::cout << sexp << alfull << std::endl << std::endl;
		std::cout << postfix << alfull << std::endl << std::endl;
		std::cout << "evals to "<< allele_t::eval(alfull.begin()) << std::endl;
		algrow.clear();
		allele_t::build(back_inserter(algrow), CREATE_GROW, i);
		std::cout << "\nGROW:\n\n";
		std::cout << sexp << algrow << std::endl << std::endl;
		std::cout << postfix << algrow << std::endl << std::endl;
		std::cout << "evals to "<< allele_t::eval(algrow.begin()) << std::endl;
		Chromosome<double> jr,sis;
		breed(alfull,algrow,std::back_inserter(jr),std::back_inserter(sis),.5);
		std::cout << sexp << "OFFSPRING: " << jr << std::endl << sis << std::endl;
	}
}
