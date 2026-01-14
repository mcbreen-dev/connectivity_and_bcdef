#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace fs {

struct Plot3DZone {
    int idx = 0; // 1-based
    std::string name;
    std::array<long long, 3> vtxSize{};
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    long long ni() const { return vtxSize[0]; }
    long long nj() const { return vtxSize[1]; }
    long long nk() const { return vtxSize[2]; }
};

std::vector<Plot3DZone> read_plot3d(const std::string& mesh_path);

} // namespace fs
