#ifndef CENO_FUNC_H
#define CENO_FUNC_H
#include <functional>
#include <memory>
#include <string>
#include <vector>

#define CENOTYPES(T) \
	using func_t = ClosedFunction<T>; \
	using func_ptr_t = std::unique_ptr<func_t>; \
	using funcs_t = std::vector<func_ptr_t>; \
	using term_t = T; \
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
	static const size_t arity = N;
};
template <class T, class... Ts> struct FuncBuilder<0, T, Ts...> {
	using type = std::function<T(Ts...)>;
};

template <class T, class... Ts> constexpr size_t arity(const std::function<T(Ts...)>&) { return sizeof...(Ts);}

template <class T> struct ClosedFunction {
	CENOTYPES(T);
	virtual size_t arity() const = 0;
	virtual term_t eval(const terms_t&) const = 0;
};

template <size_t N, class T> struct Function : public ClosedFunction<T> {
	CENOTYPES(T);
	using f_t = typename FuncBuilder<N, T>::type;
	f_t f;
	Function(const f_t& ft = f_t {}) : f(ft) {}
	constexpr size_t arity() const { return N;}
	term_t eval(const terms_t& t) const { return typename BinderBuilder<N>::binder()(f_t(), t)(); }
};

template <class T> using unary_func      = Function<1, T>;
template <class T> using binary_func     = Function<2, T>;
template <class T> using ternary_func    = Function<3, T>;
template <class T> using quaternary_func = Function<4, T>;
template <class T> using quinary_func    = Function<5, T>;

template <class T> struct Nucleotides {
	CENOTYPES(T);

	funcs_t funcs;
	terms_t terms;
	names_t func_names, term_names;
};

}// namespace ceno

#endif//CENO_FUNC_H
