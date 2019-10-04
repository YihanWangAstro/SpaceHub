//
// Created by Yihan on 4/9/19.
#include "../../src/particle-system/chain.hpp"
#include "../../src/type-class.hpp"
#include "gtest/gtest.h"

namespace UnitTest {
template <typename Coord, typename IdxArray>
void test_chain_create(Coord const &pos, IdxArray const &expected_idx) {
  IdxArray idx;
  space::Chain::calc_chain_index(pos, idx);

  IdxArray expected{4, 5, 0, 1, 2, 3};
  EXPECT_EQ(idx, expected_idx);
}

template <typename Coord, typename IdxArray>
void test_update_chain(Coord &pos, IdxArray const &idx, IdxArray const &new_idx, Coord const &expected) {
  space::Chain::update_chain(pos, idx, new_idx);
  EXPECT_EQ(pos, expected);
}
}  // namespace UnitTest
//

TEST(Chain, CreateIdx) {
  using type_sys = space::Types<double>;
  using Coord = typename type_sys::Coord;
  using IdxArray = typename type_sys::IdxArray;

  IdxArray idx{4, 5, 0, 1, 2, 3};
  Coord pos;
  pos.emplace_back(0, 0, 0);
  pos.emplace_back(1, 1, 0);
  pos.emplace_back(2, 3, 0);
  pos.emplace_back(-1, 3, 0);
  pos.emplace_back(5, -1, 0);
  pos.emplace_back(-1, -4, 0);

  UnitTest::test_chain_create(pos, idx);

  pos.clear();
  pos.emplace_back(0, 0, 0);
  pos.emplace_back(0, 1, 1);
  pos.emplace_back(0, 2, 3);
  pos.emplace_back(0, -1, 3);
  pos.emplace_back(0, 5, -1);
  pos.emplace_back(0, -1, -4);

  UnitTest::test_chain_create(pos, idx);
}