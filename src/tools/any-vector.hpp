//
// Created by yihan on 6/5/19.
//

#ifndef SPACEHUB_ANY_VECTOR_HPP
#define SPACEHUB_ANY_VECTOR_HPP

#include <any>
#include <iostream>
#include <vector>
/*---------------------------------------------------------------------------*\
        Class any Declaration
\*---------------------------------------------------------------------------*/
class AnyVector {
 public:
  AnyVector() = default;
  AnyVector(AnyVector const&) = default;
  AnyVector(AnyVector&&) = default;

  AnyVector& operator=(AnyVector const&);
  AnyVector& operator=(AnyVector&&);

  bool empty();
  size_t size();
  size_t capacity();
  void clear();
  void reserve(size_t new_cap);
  template <typename T>
  void emplace_back(T&& elem);

 private:
  std::vector<std::any> vec_;
};
/*---------------------------------------------------------------------------*\
        Class any Implementation
\*---------------------------------------------------------------------------*/
bool AnyVector::empty() { return vec_.empty(); }

size_t AnyVector::size() { return vec_.size(); }

void AnyVector::clear() { vec_.clear(); }

size_t AnyVector::capacity() { return vec_.capacity(); }

void AnyVector::reserve(size_t new_cap) { vec_.reserve(new_cap); }

template <typename T>
void AnyVector::emplace_back(T&& elem) {
  vec_.emplace_back(std::make_any(elem));
}
/*---------------------------------------------------------------------------*\
        Help function
\*---------------------------------------------------------------------------*/

std::ostream& operator<<(std::ostream& os, AnyVector const& vec) {
  for (auto const& elem : vec) {
    os << std::any_cast<decltype(elem)>(elem) << ' ';  // wrong....need to think about it.
  }
}
/*---------------------------------------------------------------------------*\
        Help macros
\*---------------------------------------------------------------------------*/
#endif  // SPACEHUB_ANY_VECTOR_HPP
