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

template <class T, class U = unsigned char> struct Allele {
	CENOTYPES(T);
	using chromosome_t = std::vector<Allele>;

	static const U TERMINAL = 1 << (std::numeric_limits<U>::digits - 1);
	static const U FUNCTION = 0;
	static Nucleotides<T> *ntides;

	U component;
	terms_t terms;
	Allele(): component(~0) {
		static_assert(std::numeric_limits<U>::is_integer && !std::numeric_limits<U>::is_signed,
			      "Index type must be unsigned");
	}
	Allele(const Allele& s): component(s.component) {}
	Allele(Allele && s): component(std::move(s.component)) {}
	Allele(const U& u): component(u) {}
	Allele& operator=(const Allele& s) {component = s.component;}
	Allele& operator=(Allele && s) {component = std::move(s.component);}
	Allele& operator=(const U& u) {component = u;}
	Allele& operator=(U && u) {component = std::move(u);}
	bool is_term() const { return component & TERMINAL;}
	bool is_func() const { return !is_term();}
	U index() const { return component & ~TERMINAL;}

	static term_t eval(typename chromosome_t::iterator iter) {
		std::stack<terms_t> argstack;
		std::stack<Allele> funcstack;
		for (; true; ++iter)
			if (iter->is_term()) {
				if (argstack.empty())
					return ntides->terms[iter->index()];
				else {
					argstack.top().push_back(ntides->terms[iter->index()]);
					while (!argstack.empty() &&
							argstack.top().size() == ntides->funcs[funcstack.top().index()]->arity()) {
						argstack.pop();
						funcstack.pop();
						if (argstack.empty())
							return ntides->funcs[funcstack.top().index()]->eval(argstack.top());
						argstack.top().push_back(ntides->funcs[funcstack.top().index()]->eval(argstack.top()));
					}
				}
			} else {
				funcstack.push(*iter);
				argstack.push(terms_t());
			}
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
		if (p & CreationPolicy::TERMS)
			if (r < ntides->terms.size()) {
				*iter = (U)(TERMINAL | r);
				++iter;
				return;
			} else
				r -= ntides->terms.size();
		*iter = (U)(FUNCTION | r);
		for (size_t i = 0; i < ntides->funcs[r]->arity(); ++i)
			build(++iter, policy, max_depth, depth + 1);
	}
};

template <class T, class U = unsigned char> using Chromosome = std::vector<Allele<T, U>>;

template <class charT, class traits, class T> inline std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& s, const ceno::Allele<T>& al) {
	if (al.is_func())
		s << (al.index() < al.ntides->func_names.size() ? al.ntides->func_names[al.index()] : "[undefined]");
	else if (al.index() < al.ntides->term_names.size())
		s << al.ntides->term_names[al.index()];
	else
		s << al.ntides->terms[al.index()];
	return s;
}

template <class charT, class traits, class T, class U> inline std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& s, const std::vector<ceno::Allele<T, U>>& als) {
	std::stack<size_t> argstack;
	std::stack<Allele<T,U>> funcstack;
	bool postfix = s.flags() & std::ios_base::right;
	for (const auto al : als) {
		if (al.is_func()) {
			argstack.push(al.ntides->funcs[al.index()]->arity());
			if (!postfix && !argstack.empty()) //no args on first pass
				s << ' ';
			if (postfix)
				funcstack.push(al);
			else
				s << '(' << al;
		}
		if (al.is_term()) {
			s << ' ' << al;
			while (!argstack.empty() && --argstack.top() == 0) {
				argstack.pop();
				if (postfix) {
					s << ' ' << funcstack.top();
					funcstack.pop();
				}else
					s << ')';
			}
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
