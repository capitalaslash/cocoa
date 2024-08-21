#pragma once

struct MeshCoupling
{
  virtual void init(
      std::string_view name,
      uint const dim,
      std::vector<double> const & coords,
      std::vector<uint> conn,
      std::vector<uint> offsets) {};
};
