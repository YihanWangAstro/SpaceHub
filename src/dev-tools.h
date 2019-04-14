#ifndef DEVTOOLS_h
#define DEVTOOLS_h

#include <iostream>
#include <tuple>
namespace space {

    template<typename... Args>
    void print(std::ostream &os, Args &&... args) {
        (os << ...<< std::forward<Args>(args));
    }

    template<typename... Args>
    void input(std::istream &is, Args &&... args) {
        (is >> ... >> std::forward<Args>(args));
    }

    template<typename... Args>
    void display(std::ostream &out, Args &&... args) {
        (..., (out << std::forward<Args>(args) << ' '));
    }

    template<class Tup, size_t... I>
    void print_tuple(std::ostream &out, const Tup &&tup, std::index_sequence<I...>) {
        (..., (out << std::get<I>(tup)));
    }

    template<class Tup, size_t... I>
    void input_tuple(std::istream &in, const Tup &&tup, std::index_sequence<I...>) {
        (..., (in >> std::get<I>(tup)));
    }

    template<typename ...Args>
    std::ostream &operator<<(std::ostream &out, std::tuple<Args...> &&tup) {
        space::print_tuple(out, std::forward<decltype(tup)>(tup), std::make_index_sequence<sizeof...(Args)>());
        return out;
    }

    template<typename ...Args>
    std::istream &operator>>(std::istream &in, std::tuple<Args...> &&tup) {
        space::input_tuple(in, std::forward<decltype(tup)>(tup), std::make_index_sequence<sizeof...(Args)>());
        return in;
    }

    template <typename STL, typename ...Args>
    void emplace_back(STL& container, Args&&...args){
        (container.emplace_back(std::forward<Args>(args)),...);
    }

    template <typename STL, typename ...Args>
    void push_back(STL& container, Args&&...args){
        (container.emplace_back(std::forward<Args>(args)),...);
    }

    template<typename... Args>
    void resize_all(size_t new_sz, Args &&... args) {
        (..., (args.resize(new_sz)));
    }

    template<typename... Args>
    void reserve_all(size_t new_cap, Args &&... args) {
        (..., (args.reserve(new_cap)));
    }

    class Empty{};

    template<typename T>
    struct get_value_type {
    private:
        /*If U has member::value_type, getValueType<T>(0) will match this function. See details on SFINAE. */
        template<typename U>
        static typename U::value_type check(typename U::value_type);

        /*If U doesn't have member::value_type, getValueType<T>(0) will match this function. See details on SFINAE. */
        template<typename U>
        static U check(U);

    public:
        using type = decltype(check<T>(0));
    };

#define MACRO_CAT(A, B) MACRO_CAT_I(A, B)
#define MACRO_CAT_I(A, B) MACRO_CAT_II(~, A ## B)
#define MACRO_CAT_II(P, REST) REST
#define UNIQ(BASE) MACRO_CAT(BASE, __LINE__)

#define SPACEHUB_ABORT(...) {                                                                                          \
    space::print(std::cout, __FILE__, ": Line :",  __LINE__ , "\r\n");                                                 \
    space::print(std::cout, __VA_ARGS__);                                                                             \
    exit(0);                                                                                                           \
}

#define PACK(...) std::forward_as_tuple(__VA_ARGS__)

#ifdef DEBUG
#define DEBUG_MODE(BLOCK) BLOCK
#else
#define DEBUG_MODE(BLOCK)
#endif

/** @brief Macros used to output debuf info.  */
#ifdef DEBUG
#define DEBUG_MSG(EXPR,...) (EXPR ? SpaceH::print(std::cout,  __VA_ARGS__ ) : void(0) )
#else
#define DEBUG_MSG(EXPR, ...)
#endif

#ifdef DEBUG
#define DEBUG_MODE_ASSERT(EXPR,MSG) ((EXPR) ? void(0) : (SPACEHUB_ABORT(MSG)))
#else
#define DEBUG_MODE_ASSERT(EXPR, MSG)
#endif

#define SPACEHUB_USING_TYPE_SYSTEM_OF(CLASS)                                                                           \
    template<typename ..._T_>                                                                                          \
    using Container   = typename CLASS::template Container<_T_...>;                                                    \
                                                                                                                       \
    using Scalar      = typename CLASS::Scalar;                                                                        \
    using ScalarArray = typename CLASS::ScalarArray;                                                                   \
    using IdxArray    = typename CLASS::IdxArray;                                                                      \
    using IntArray    = typename CLASS::IntArray;                                                                      \
    using Vector      = typename CLASS::Vector;                                                                        \
    using VectorArray = typename CLASS::VectorArray;                                                                   \
    using Coord       = typename CLASS::Coord

#define DECLARE_CRTP_ACCESSOR(DERIVED, TYPE, NAME)                                                                     \
                                                                                                                       \
inline TYPE & NAME () {                                                                                                \
    return static_cast<Derived*>(this)->impl_##NAME();                                                                 \
};                                                                                                                     \
inline TYPE const & NAME () const {                                                                                    \
    return static_cast<Derived const*>(this)->impl_##NAME();                                                           \
};

#define DECLARE_READ_ACCESSOR(NAME, TYPE, DERIVED)                                                                     \
                                                                                                                       \
inline TYPE const & NAME () const {                                                                                    \
    return static_cast<Derived*>(this)->impl_##NAME();                                                                 \
};

#define SPACEHUB_STD_ACCESSOR(TYPE, NAME, MEMBER)                                                                      \
                                                                                                                       \
inline TYPE & NAME () {                                                                                                \
    return MEMBER;                                                                                                     \
};                                                                                                                     \
inline TYPE const & NAME () const {                                                                                    \
    return MEMBER;                                                                                                     \
};

#define SPACEHUB_READ_ACCESSOR(TYPE, NAME, MEMBER)                                                                     \
                                                                                                                       \
inline TYPE const & NAME () const {                                                                                    \
    return MEMBER;                                                                                                     \
};

#define CREATE_METHOD_CHECK(NAME)                                                                                      \
                                                                                                                       \
    template<typename T, typename... Args>                                                                             \
    struct has_method_##NAME                                                                                           \
    {                                                                                                                  \
        template<typename U>                                                                                           \
        constexpr static auto check(const void*)                                                                       \
        ->decltype(std::declval<U>().NAME(std::declval<Args>()...), std::true_type());                                 \
                                                                                                                       \
        template<typename U>                                                                                           \
        constexpr static std::false_type check(...);                                                                   \
                                                                                                                       \
        static constexpr bool value = decltype(check<T>(nullptr))::value;                                              \
    };

#define HAS_METHOD(CLASS, METHOD, ...) has_method_##METHOD<CLASS, ##__VA_ARGS__>::value

/** @brief Macros used to check if a class has a member.*/
#define CREATE_MEMBER_CHECK(MEMBER)                                                                                    \
                                                                                                                       \
    template<typename T, typename V = bool>                                                                            \
    struct has_ ## MEMBER : std::false_type { };                                                                       \
                                                                                                                       \
    template<typename T>                                                                                               \
    struct has_ ## MEMBER                                                                                              \
    <                                                                                                                  \
        T,                                                                                                             \
        typename std::enable_if                                                                                        \
        <                                                                                                              \
            !std::is_same<decltype(std::declval<T>().MEMBER), void>::value, bool                                       \
        >::type                                                                                                        \
    >: std::true_type { };

#define HAS_MEMBER(C, member) has_ ## member<C>::value

/** @brief Macros used to check if a class has a static member.  */
#define CREATE_STATIC_MEMBER_CHECK(MEMBER)                                                                             \
                                                                                                                       \
    template<typename T, typename V = bool>                                                                            \
    struct has_static ## MEMBER : std::false_type { };                                                                 \
                                                                                                                       \
    template<typename T>                                                                                               \
    struct has_static ## MEMBER                                                                                        \
    <                                                                                                                  \
        T,                                                                                                             \
        typename std::enable_if                                                                                        \
        <                                                                                                              \
            !std::is_same<decltype(T::MEMBER), void>::value, bool                                                      \
        >::type                                                                                                        \
    >: std::true_type { };

#define HAS_STATIC_MEMBER(C, member) has_static ## member<C>::value

/** @brief Macros used to static_assert if a class has a specific method. */
#define CHECK_METHOD(CLASS, METHOD, ...)                                                                               \
                                                                                                                       \
            static_assert(has_method_##METHOD<CLASS, ##__VA_ARGS__>::value,                                            \
            "Template argument '" # CLASS  "' must have method '"  # METHOD  "("  #__VA_ARGS__  ")'. ");

/** @brief Macros used to static_assert if a class has a member. */
#define CHECK_MEMBER(CLASS, MB)                                                                                        \
                                                                                                                       \
            static_assert(has_ ## MB<C>::value,                                                                        \
            "Template argument '" # CLASS  "' must have member '"  # MB  "'. ");

/** @brief Macros used to static_assert if a class has a static member.*/
#define CHECK_STATIC_MEMBER(CLASS, MB)                                                                                 \
                                                                                                                       \
            static_assert(has_static ## MB<CLASS>::value,                                                              \
            "Template argument '" # CLASS  "' must have static member '"  # MB  "'. ");

/** @brief Macros used to static_assert if two class have the same base type set*/
#define CHECK_TYPE(T1, T2)                                                                                             \
            static_assert(std::is_same< typename T1::Types, typename T2::Types>::value,                                \
            "Template argument '" #T1 "' and '" #T2 "' must have the same type of the type member(SpaceH::ProtoType<...>)");


#define CHECK_POD(DATA)                                                                                                \
            static_assert(std::is_trivial<DATA>::value, "Template arg '" #DATA "' must be a POD type!");

    template<typename T, typename... Args>
    struct is_indexable {
        template<typename U>
        constexpr static auto check(const void *)
        -> decltype(std::declval<U>().operator[](std::declval<Args>()...), std::true_type());

        template<typename U>
        constexpr static std::false_type check(...);

        static constexpr bool value = decltype(check<T>(nullptr))::value;
    };

    template<typename T, typename _ = void>
    struct is_container : std::false_type {};

    template<typename... Ts>
    struct is_container_helper {};

    template<typename T>
    struct is_container<
            T,
            std::conditional_t<
                    false,
                    is_container_helper<
                            typename T::value_type,
                            typename T::size_type,
                            //typename T::allocator_type,
                            typename T::iterator,
                            typename T::const_iterator,
                            decltype(std::declval<T>().size()),
                            decltype(std::declval<T>().begin()),
                            decltype(std::declval<T>().end())
                            //decltype(std::declval<T>().cbegin()),
                            //decltype(std::declval<T>().cend())
                    >,
                    void
            >
    > : public std::true_type {};

    template <typename T>
    constexpr bool is_container_v = is_container<T>::value;

//#define INDEXABLE(TYPE) ((is_indexable<TYPE, size_t>::value)||(is_indexable<TYPE, int>::value))

#define IS_BASE_OF(BASE, DERIVED) (std::is_base_of<BASE,DERIVED>::value)
#define TYPE_OF_SELF std::remove_reference<decltype(*this)>::type
#define REF_TYPE_OF_SELF decltype(*this)
}//end namespace SpaceH

#endif
