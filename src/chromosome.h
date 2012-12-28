#ifndef CENO_allele_H
#define CENO_allele_H
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <stack>
#include <vector>

#include "rand.h"
#include "nucleotide.h"

namespace ceno {

template <class T, class U = unsigned char> class Allele {
	CENOTYPES(T);
	using genome_t = std::vector<Allele>;
	using locus_t = typename genome_t::const_iterator;
	using chromosome_t = locus_t[2];

	static const U TERMINAL = 1 << (std::numeric_limits<U>::digits - 1);
	static const U FUNCTION = 0;
	U component = ~0;
	public:
	static Nucleotides<T> *ntides;
	Allele() {
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

	static term_t eval(typename genome_t::iterator iter) {
		std::stack<terms_t> argstack;
		std::stack<Allele> funcstack;
		for (; true; ++iter)
			if (iter->is_term()) {
				if (argstack.empty())
					return ntides->terms[iter->index()];
				else {
					argstack.top().push_back(ntides->terms[iter->index()]);
					for (auto &a : argstack.top()) {
					}
					while ( argstack.top().size() == ntides->funcs[funcstack.top().index()]->arity()) {
						term_t t = ntides->funcs[funcstack.top().index()]->eval(argstack.top());
						funcstack.pop();
						argstack.pop();
						if (argstack.empty())
							return t;
						argstack.top().push_back(t);
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

	static locus_t random(const locus_t& begin, const locus_t& end, bool wantfunc) {
		size_t r = bounded_rand(end-begin);
		if ((begin+r)->is_func() ==  wantfunc)
			return begin + r;
		else {
			auto iter = begin + r;
			//every terminal has to have a preceding function, every function a following terminal
			for (;iter->is_func() != wantfunc; iter += wantfunc ? -1 : 1) /*empty*/;
			return iter;
		}
	}

	static locus_t extract(locus_t& l) {
		std::stack<size_t> argstack;
		if (l->is_term())
			return ++l;
		for (; true; ++l ){
			if (l->is_func())
				argstack.push(l->ntides->funcs[l->index()]->arity());
			if (l->is_term())
				while (!argstack.empty() && --argstack.top() == 0)
					argstack.pop();
			if (argstack.empty())
				return ++l;
		}
	}

};

template <class T, class U> using Genome = std::vector<Allele<T, U>>;
template <class T, class U> using Locus = typename Genome<T, U>::const_iterator;

template <class T, class U>  class Chromosome {
	using allele_t = Allele<T,U>;
	using locus_t = Locus<T,U>;
	using genome_t = Genome<T,U>;
    locus_t b,e;
	public:
	Chromosome(const locus_t & l):b(l){ e = allele_t::extract(b); }
	Chromosome(const genome_t& g):b(std::begin(g)), e(std::end(g)) {}
	Chromosome& operator=(const genome_t& g) { b = std::begin(g); e = std::end(g); }
	const locus_t& begin() const { return b; }
	const locus_t& end() const { return e; }
	size_t size() const { return e - b; }
	static Chromosome random(const Chromosome &c, float functionweight){
		return Chromosome(allele_t::random(std::begin(c),std::end(c),bounded_rand(10000) < functionweight * 10000));
	}
};

template <class T, class U, class Out>
void breed(const Genome<T,U> &xx, const Genome<T,U> &xy, Out jr, Out sis, float functionweight = .8){
	auto fromdad = Chromosome<T,U>::random(xy,functionweight);
	auto frommom = Chromosome<T,U>::random(xx,functionweight);
	std::copy(xy.begin(),fromdad.begin(),jr);
	std::copy(frommom.begin(),frommom.end(),jr);
	std::copy(fromdad.end(),xy.end(),jr);
	std::copy(xx.begin(),frommom.begin(),sis);
	std::copy(fromdad.begin(),fromdad.end(),sis);
	std::copy(frommom.end(),xx.end(),sis);
}

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
