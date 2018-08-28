#ifndef DEVTOOLS_h
#define DEVTOOLS_h

#include<iostream>

namespace SpaceH {
    /**
     *
     * @tparam Arg
     * @tparam Args
     * @param out
     * @param arg
     * @param args
     */
    template<typename Arg, typename... Args>
    void print(std::ostream &out, Arg &&arg, Args &&... args) {
        out << std::forward<Arg>(arg);
        using expander = int[];
        (void) expander{0, (void(out << ' ' << std::forward<Args>(args)), 0)...};
        out << '\n';
    }

/**
 *
 * @param msg
 * @param file
 * @param line
 */
    void errMsg(const char *msg, const char *file, size_t line) {
        std::cout << "An error occurred:" << '\n';
        std::cout << "  Message >>>" << msg << '\n';
        std::cout << "  File    >>>" << file << '\n';
        std::cout << "  Line    >>>" << line << std::endl;
        exit(0);
    }


/** @brief print an array. Used for debug*/
    template<typename T>
    void printArray(T &var) {
        for (const auto & v : var)
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

/** @brief Macros used to output debuf info.  */
#ifdef DEBUG
#define DEBUG_MSG(cond,...) ( cond ? SpaceH::print(std::cout,  __VA_ARGS__ ) : void(0) )
#else
#define DEBUG_MSG(cond, ...)
#endif

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
            static_assert(std::is_same< typename T1::type, typename T2::type>::value,               \
            "Template argument '" #T1 "' and '" #T2 "' must have the same type of the type member(SpaceH::ProtoType<...>)");


#define CHECK_POD(DATA)                                                                             \
            static_assert(std::is_trivial<DATA>::value, "Template arg '" #DATA "' must be a POD type!");
}//end namespace SpaceH

#endif
