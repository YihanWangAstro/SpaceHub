#ifndef DEVTOOLS_h
#define DEVTOOLS_h

#include<iostream>

namespace SpaceH {

    template<typename Arg, typename... Args>
    void print(std::ostream &out, Arg &&arg, Args &&... args) {
        out << std::forward<Arg>(arg);
        using expander = int[];
        (void) expander{0, (void(out << ' ' << std::forward<Args>(args)), 0)...};
        out << '\n';
    }

#define SPACEHUB_ERR_MSG(...) {                                 \
	SpaceH::print(std::cout, __FILE__, ": Line :",  __LINE__ ); \
	SpaceH::print(std::cout, __VA_ARGS__ );                     \
	exit(0);                                                    \
}


/** @brief print an array. Used for debug*/
    template<typename T>
    void printArray(T &var) {
        for (const auto &v : var)
            std::cout << v << '\n';

        std::cout << '\n';
    }

/**
 *
 * @tparam T
 */
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

    template<typename T>
    struct get_types {
    private:
        template<typename U, typename V = bool>
        struct has_types : std::false_type {};

        template<typename U>
        struct has_types<U, typename std::enable_if<!std::is_same<typename U::Types, void>::value, bool>::type> : std::true_type {};

        /*If U has member::Types, check<T>() will match this function. See details on SFINAE. */
        template<typename X>
        static typename std::enable_if<has_types<X>::value, typename X::Types>::type check();

        /*If U doesn't have member::Types, check<T>() will match this function. See details on SFINAE. */
        template<typename X>
        static typename std::enable_if<!has_types<X>::value, X>::type check();

    public:
        using Types = decltype(check<T>());
    };

/** @brief Macros used to output debuf info.  */
#ifdef DEBUG
#define DEBUG_MSG(cond,...) ( cond ? SpaceH::print(std::cout,  __VA_ARGS__ ) : void(0) )
#else
#define DEBUG_MSG(cond, ...)
#endif

#define SPACEHUB_USING_TYPE_SYSTEM_OF(CLASS)                                                              \
    using Types        = typename SpaceH::get_types<CLASS>::Types;                                        \
    using Scalar       = typename CLASS::Scalar;                                                          \
    using Vector       = typename CLASS::Vector;                                                          \
    using VectorArray  = typename CLASS::VectorArray;                                                     \
    using ScalarArray  = typename CLASS::ScalarArray;                                                     \
    using IndexArray   = typename CLASS::IndexArray;                                                      \
    using ScalarBuffer = typename CLASS::ScalarBuffer;                                                    \
    using SizeArray    = typename CLASS::SizeArray;                                                       \
    template<typename T, size_t S>                                                                        \
    using Container    = typename CLASS::template Container<T, S>;

/** @brief Standard read interfaces for private data scalar in SpaceHub project*/
#define SPACEHUB_READ_INTERFACES_FOR_SCALAR(NAME, TYPE, MEMBER)                                          \
                                                                                                         \
inline const TYPE & NAME () const {                                                                      \
    return MEMBER;                                                                                       \
};

/** @brief Standard write interfaces for private data scalar in SpaceHub project*/
#define SPACEHUB_WRITE_INTERFACES_FOR_SCALAR(NAME, TYPE, MEMBER)                                         \
inline void set_##NAME (const TYPE& scalar) {                                                            \
    MEMBER = scalar;                                                                                     \
};                                                                                                       \
inline void swap_##NAME (TYPE& scalar) {                                                                 \
    std::swap(MEMBER, scalar);                                                                           \
}

/** @brief Standard scalar read interfaces adapter in SpaceHub project*/
#define SPACEHUB_READ_INTERFACES_ADAPTER_FOR_SCALAR(NEWNAME, TYPE, MEMBER, NAME)                         \
                                                                                                         \
inline const TYPE & NEWNAME () const {                                                                   \
    return MEMBER.NAME();                                                                                \
};


/** @brief Standard scalar read interfaces adapter in SpaceHub project*/
#define SPACEHUB_READ_INTERFACES_ADAPTER_FOR_BUILDIN(NEWNAME, TYPE, MEMBER, NAME)                        \
                                                                                                         \
inline const TYPE & NEWNAME () const {                                                                   \
    return MEMBER.NAME();                                                                                \
};

/** @brief Standard scalar write interfaces adapter in SpaceHub project*/
#define SPACEHUB_WRITE_INTERFACES_ADAPTER_FOR_SCALAR(NEWNAME, TYPE, MEMBER, NAME)                        \
inline void set_##NEWNAME (const TYPE& scalar) {                                                         \
    MEMBER.set_##NAME(scalar);                                                                           \
};                                                                                                       \
inline void swap_##NEWNAME (TYPE& scalar) {                                                              \
    MEMBER.swap_##NAME(scalar);                                                                          \
};

/** @brief Standard read interfaces for private data array in SpaceHub project*/
#define SPACEHUB_READ_INTERFACES_FOR_ARRAY(NAME, TYPE, MEMBER)                                           \
                                                                                                         \
inline const TYPE & NAME () const {                                                                      \
    return MEMBER;                                                                                       \
};                                                                                                       \
inline const typename TYPE::value_type & NAME (size_t i) const {                                         \
    return MEMBER[i];                                                                                    \
};

/** @brief Standard write interfaces for private data array in SpaceHub project*/
#define SPACEHUB_WRITE_INTERFACES_FOR_ARRAY(NAME, TYPE, MEMBER)                                          \
inline void set_##NAME (const TYPE& array) {                                                             \
    MEMBER = array;                                                                                      \
};                                                                                                       \
inline void set_##NAME (size_t i, const typename TYPE::value_type &value) {                              \
    MEMBER[i] = value;                                                                                   \
};                                                                                                       \
inline void swap_##NAME (TYPE& array) {                                                                  \
    std::swap(array, MEMBER);                                                                            \
};

/** @brief Standard read array interfaces adapter in SpaceHub project*/
#define SPACEHUB_READ_INTERFACES_ADAPTER_FOR_ARRAY(NEWNAME, TYPE, MEMBER, NAME)                          \
                                                                                                         \
inline const TYPE & NEWNAME () const {                                                                   \
    return MEMBER.NAME();                                                                                \
};                                                                                                       \
inline const typename TYPE::value_type & NEWNAME (size_t i) const {                                      \
    return MEMBER.NAME(i);                                                                               \
};

/** @brief Standard write array interfaces adapter in SpaceHub project*/
#define SPACEHUB_WRITE_INTERFACES_ADAPTER_FOR_ARRAY(NEWNAME, TYPE, MEMBER, NAME)                         \
inline void set_##NEWNAME (const TYPE& array) {                                                          \
    MEMBER.set_##NAME(array);                                                                            \
};                                                                                                       \
inline void set_##NEWNAME (size_t i, const typename TYPE::value_type &value) {                           \
    MEMBER.set_##NAME(i, value);                                                                         \
};                                                                                                       \
inline void swap_##NEWNAME (TYPE& array) {                                                               \
    MEMBER.swap_##NAME(array);                                                                           \
};

/** @brief Macros used to check if a class has a specific method.  */
#define CREATE_METHOD_CHECK(NAME)                                                                   \
                                                                                                    \
    template<typename T, typename... Args>                                                          \
    struct has_method_##NAME                                                                        \
    {                                                                                               \
        template<typename U>                                                                        \
        constexpr static auto check(const void*)                                                    \
        ->decltype(std::declval<U>().NAME(std::declval<Args>()...), std::true_type());              \
                                                                                                    \
        template<typename U>                                                                        \
        constexpr static std::false_type check(...);                                                \
                                                                                                    \
        static constexpr bool value = decltype(check<T>(nullptr))::value;                           \
    };

#define HAS_METHOD(CLASS, METHOD, ...) has_method_##METHOD<CLASS, __VA_ARGS__>::value

/** @brief Macros used to check if a class has a member.  */
#define CREATE_MEMBER_CHECK(MEMBER)                                                                 \
                                                                                                    \
    template<typename T, typename V = bool>                                                         \
    struct has_ ## MEMBER : std::false_type { };                                                    \
                                                                                                    \
    template<typename T>                                                                            \
    struct has_ ## MEMBER                                                                           \
    <                                                                                               \
        T,                                                                                          \
        typename std::enable_if                                                                     \
        <                                                                                           \
            !std::is_same<decltype(std::declval<T>().MEMBER), void>::value, bool                    \
        >::type                                                                                     \
    >: std::true_type { };

#define HAS_MEMBER(C, member) has_ ## member<C>::value

/** @brief Macros used to check if a class has a static member.  */
#define CREATE_STATIC_MEMBER_CHECK(MEMBER)                                                          \
                                                                                                    \
    template<typename T, typename V = bool>                                                         \
    struct has_static ## MEMBER : std::false_type { };                                              \
                                                                                                    \
    template<typename T>                                                                            \
    struct has_static ## MEMBER                                                                     \
    <                                                                                               \
        T,                                                                                          \
        typename std::enable_if                                                                     \
        <                                                                                           \
            !std::is_same<decltype(T::MEMBER), void>::value, bool                                   \
        >::type                                                                                     \
    >: std::true_type { };

#define HAS_STATIC_MEMBER(C, member) has_static ## member<C>::value


/** @brief Macros used to static_assert if a class has a specific method. */
#define CHECK_METHOD(CLASS, METHOD, ...)                                                            \
                                                                                                    \
            static_assert(has_method_##METHOD<CLASS, ##__VA_ARGS__>::value,                         \
            "Template argument '" # CLASS  "' must have method '"  # METHOD  "("  #__VA_ARGS__  ")'. ");

/** @brief Macros used to static_assert if a class has a member. */
#define CHECK_MEMBER(CLASS, MB)                                                                     \
                                                                                                    \
            static_assert(has_ ## MB<C>::value,                                                     \
            "Template argument '" # CLASS  "' must have member '"  # MB  "'. ");

/** @brief Macros used to static_assert if a class has a static member.*/
#define CHECK_STATIC_MEMBER(CLASS, MB)                                                              \
                                                                                                    \
            static_assert(has_static ## MB<CLASS>::value,                                           \
            "Template argument '" # CLASS  "' must have static member '"  # MB  "'. ");

/** @brief Macros used to static_assert if two class have the same base type set*/
#define CHECK_TYPE(T1, T2)                                                                          \
            static_assert(std::is_same< typename T1::Types, typename T2::Types>::value,             \
            "Template argument '" #T1 "' and '" #T2 "' must have the same type of the type member(SpaceH::ProtoType<...>)");


#define CHECK_POD(DATA)                                                                             \
            static_assert(std::is_trivial<DATA>::value, "Template arg '" #DATA "' must be a POD type!");
}//end namespace SpaceH

#endif
