#ifndef CENO_allele_H
#define CENO_allele_H
#include <cassert>
#include <iostream>
#include <limits>
#include <stack>
#include <vector>

#include "rand.h"
#include "nucleotide.h"

namespace ceno {

template <class F, class T = unsigned char> struct Allele {
	CENOTYPES(F);
	using snodes_t = std::vector<Allele>;

	static const T TERMINAL = 1 << (std::numeric_limits<T>::digits - 1);
	static const T FUNCTION = 0;
	static Nucleotides<func_t> *ntides;

	T component;
	terms_t terms;
	Allele(): component(~0) {}
	Allele(const Allele& s): component(s.component) {}
	Allele(Allele && s): component(std::move(s.component)) {}
	Allele(const T& t): component(t) {}
	Allele& operator=(const Allele& s) {component = s.component;}
	Allele& operator=(Allele && s) {component = std::move(s.component);}
	Allele& operator=(const T& t) {component = t;}
	Allele& operator=(T && t) {component = std::move(t);}
	bool is_term() const { return component & TERMINAL;}
	bool is_func() const { return !is_term();}
	T index() const { return component & ~TERMINAL;}

	term_t eval(typename snodes_t::iterator iter) const {
		if (is_term())
			return ntides->terminals[index()];
		std::array<term_t, arity(func_t())> args;
		for (int i = 0; i < args.size(); ++i)
			args[i] = eval(++iter);
		typename BinderBuilder<args.size()>::binder b;
		auto bound = b(ntides->functions[index()], args);
		return bound();
	}

	template <class Iter>
	static void build(Iter iter, int policy, size_t max_depth, size_t depth = 0) {
		auto p = policy;
		if (depth == 0)
			p = CreationPolicy::FUNCS;
		if (depth >= max_depth)
			p = CreationPolicy::TERMS;
		size_t bound = ((p & CreationPolicy::TERMS) ? ntides->terms.size() : 0)
			       + ((p & CreationPolicy::FUNCS) ? ntides->funcs.size() : 0);
		auto r = bounded_rand(bound);
		//std::cerr << "bound == " << bound << " r == " << r << std::endl;
		if (p & CreationPolicy::TERMS)
			if (r < ntides->terms.size()) {
				*iter = (T)(TERMINAL | r);
				//std::cerr << "term == " << (TERMINAL | r) << std::endl;
				++iter;
				return;
			} else
				r -= ntides->terms.size();
		*iter = (T)(FUNCTION | r);
		//std::cerr << "func == " << (FUNCTION | r) << std::endl;
		for (size_t i = 0; i < arity(func_t()); ++i)
			build(++iter, policy, max_depth, depth + 1);
	}
};

template <class T> using Chromosome = std::vector<Allele<T>>;

template <class charT, class traits, class F> inline std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& s, const ceno::Allele<F>& al) {
	if (al.is_func())
		s << (al.index() < al.ntides->func_names.size() ? al.ntides->func_names[al.index()] : "[undefined]");
	else if (al.index() < al.ntides->term_names.size())
		s << al.ntides->term_names[al.index()];
	else
		s << al.ntides->terms[al.index()];
	return s;
}

template <class charT, class traits, class F> inline std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& s, const std::vector<ceno::Allele<F>>& als) {
	std::stack<size_t> argstack;
for (const auto al : als) {
		if (!argstack.empty()) //no args on first pass
			s << ' ';
		if (al.is_func()) {
			s << '(';
			argstack.push(arity(F()));
		}
		s << al;
		if (al.is_term())
			while (!argstack.empty() && --argstack.top() == 0) {
				argstack.pop();
				s << ')';
			}
	}
	assert(argstack.empty());
	return s;
}

std::ios_base& prefix(std::ios_base& str);
std::ios_base& postfix(std::ios_base& str);
std::ios_base& sexp(std::ios_base& str);
std::ios_base& forth(std::ios_base& str);

}//namespace ceno

#endif//CENO_allele_H
