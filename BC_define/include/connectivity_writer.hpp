/*
  Connectivity-writer declarations
*/
#pragma once

#include "common.hpp"
#include <cstdint>

namespace fs {

class Mesh;

struct ConnPatch
{
    int      recvZone = 0;
    int      donorZone = 0;
    FaceDir  recvDir = FaceDir::IMIN;
    FaceDir  donorDir = FaceDir::IMIN;
    long long offU = 0;
    long long offV = 0;
    long long sizeU = 0;
    long long sizeV = 0;
    uint8_t  ori = 0;
    int      transform[3] = {0, 0, 0};
    PointRange recvRange{};
    PointRange donorRange{};
    bool    useExplicitRanges = false;
};

extern const int map[8][3];

void write_1to1(Mesh& mesh, const ConnPatch& cp);

} // namespace fs
