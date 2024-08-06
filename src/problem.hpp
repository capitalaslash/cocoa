#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>

struct MEDField;

struct Problem
{
  using ParamList_T = std::unordered_map<std::string, std::filesystem::path>;

  Problem() = default;
  virtual ~Problem() = default;

  virtual void setup(ParamList_T const & params) = 0;
  virtual bool run() = 0;
  virtual void advance() = 0;
  virtual void solve() = 0;
  virtual void print() = 0;
  virtual MEDField get_field(std::string_view name) = 0;
  virtual void set_field(std::string_view name, MEDField const & field) = 0;

  double time = 0.0;
  int it = 0;
};
