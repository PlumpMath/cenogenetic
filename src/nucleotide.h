#ifndef CENO_FUNC_H
#define CENO_FUNC_H
#include <functional>
#include <string>
#include <vector>

#define CENOTYPES(F) \
	using func_t = F; \
	using funcs_t = std::vector<func_t>; \
	using term_t = typename func_t::result_type; \
	using terms_t = std::vector<term_t>; \
	using names_t = std::vector<std::string>;

namespace ceno {

enum CreationPolicy {
	FUNCS		= 1 << 0,
	ARGS		= 1 << 1,
	CONSTS		= 1 << 2,
	TERMS		= CONSTS | ARGS,
	ALL			= FUNCS | TERMS
};
const int CREATE_FULL = CreationPolicy::FUNCS;
const int CREATE_GROW = CreationPolicy::ALL;

template <size_t... Indices> struct Binder {
	template <typename F, typename C> auto operator()(F f, const C& c)
	-> decltype(std::bind(f, c.at(Indices)...)) {
		return std::bind(f, c.at(Indices)...);
	}
};

template <size_t X, size_t Y = 0, size_t... Inds> struct BinderBuilder {
	using binder = typename BinderBuilder < X - 1, Y + 1, Y, Inds... >::binder;
};
template <size_t Y, size_t... Inds> struct BinderBuilder<0, Y, Inds...> {
	using binder = Binder<Inds...>;
};

template <size_t N, class T, class... Ts> struct FuncBuilder {
	using type = typename FuncBuilder < N - 1, T, T, Ts... >::type;
};
template <class T, class... Ts> struct FuncBuilder<0, T, Ts...> {
	using type = std::function<T(Ts...)>;
};

template <class T, class... Ts> constexpr size_t arity(const std::function<T(Ts...)>&) { return sizeof...(Ts);}

template <class T> using unary_func      = typename FuncBuilder<1, T>::type;
template <class T> using binary_func     = typename FuncBuilder<2, T>::type;
template <class T> using ternary_func    = typename FuncBuilder<3, T>::type;
template <class T> using quaternary_func = typename FuncBuilder<4, T>::type;
template <class T> using quinary_func    = typename FuncBuilder<5, T>::type;

template <class F> struct Nucleotides {
	CENOTYPES(F);

	funcs_t funcs;
	terms_t terms;
	names_t func_names, term_names;
};

}// namespace ceno

#endif//CENO_FUNC_H
