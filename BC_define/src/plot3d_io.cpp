/*─────────────────────────────────────────────────────────────
  File: src/plot3d_io.cpp

  Plot3D mesh reader (Fortran unformatted, multi-block).

  This file implements a minimal reader for structured Plot3D grids
  written in classic Fortran unformatted (record-based) format.

  Supported format assumptions:
    - Single file containing:
        * Record 1: number of blocks (int32)
        * Record 2: dimensions for each block (ni,nj,nk) as int32 triplets
        * Record 3+: one record per block containing coordinates:
            [x(1..N), y(1..N), z(1..N)] as double precision
    - Coordinates are stored block-by-block
    - Endianness matches the host (no byte swapping performed)

  The reader produces a vector<Plot3DZone>, which is later consumed
  by the connectivity detection pipeline (connectivity_core.cpp).

  Notes:
    - No flow variables are read here (grid-only Plot3D).
    - No blanking, iblank arrays, or solution files are handled.
    - This is intentionally strict: unexpected record sizes throw.
─────────────────────────────────────────────────────────────*/
#include "plot3d_io.hpp"

#include <cstring>
#include <fstream>
#include <stdexcept>

namespace {

/*=====================================================================
  FortranRecordReader

  Helper class for reading Fortran unformatted sequential records.

  A Fortran unformatted record has the structure:
      [uint32 record_size]
      [record_size bytes of payload]
      [uint32 record_size]

  This class:
    - Validates matching leading/trailing record markers
    - Returns raw payload bytes as std::vector<char>
    - Provides a convenience method for reading scalar records
=====================================================================*/
class FortranRecordReader {
public:
    explicit FortranRecordReader(const std::string& path)
        : input_(path, std::ios::binary) {
        if (!input_) {
            throw std::runtime_error("Failed to open input file: " + path);
        }
    }

    /*---------------------------------------------------------
      Read a single Fortran record and return its payload.
    ---------------------------------------------------------*/
    std::vector<char> readRecord() {
        std::uint32_t leading = 0;
        if (!input_.read(reinterpret_cast<char*>(&leading), sizeof(leading))) {
            throw std::runtime_error("Unable to read record marker from file.");
        }

        std::vector<char> buffer(leading);
        if (!buffer.empty() && !input_.read(buffer.data(), leading)) {
            throw std::runtime_error("Unable to read record payload from file.");
        }

        std::uint32_t trailing = 0;
        if (!input_.read(reinterpret_cast<char*>(&trailing), sizeof(trailing))) {
            throw std::runtime_error("Unable to read trailing record marker from file.");
        }

        if (leading != trailing) {
            throw std::runtime_error("Mismatched Fortran record markers detected.");
        }

        return buffer;
    }

    /*---------------------------------------------------------
      Read a scalar value stored as a Fortran record.

      The record payload size must exactly match sizeof(T).
    ---------------------------------------------------------*/
    template <typename T>
    T readScalar() {
        auto buffer = readRecord();
        if (buffer.size() != sizeof(T)) {
            throw std::runtime_error("Unexpected record length while reading scalar.");
        }
        T value{};
        std::memcpy(&value, buffer.data(), sizeof(T));
        return value;
    }

private:
    std::ifstream input_;
};

} // anonymous namespace

namespace bcdef {

/*=====================================================================
  read_plot3d

  Read a multi-block structured Plot3D grid file.

  Parameters:
    mesh_path : path to the Plot3D grid file

  Returns:
    vector<Plot3DZone>, one entry per block, populated with:
      - idx      : 1-based block index
      - name     : "block-N"
      - vtxSize  : (ni,nj,nk)
      - x,y,z    : coordinate arrays, flattened in Plot3D order

  Layout assumptions:
    - Coordinates are stored as:
          x(1..N), y(1..N), z(1..N)
      where N = ni * nj * nk
=====================================================================*/
std::vector<Plot3DZone> read_plot3d(const std::string& mesh_path) {
    FortranRecordReader reader(mesh_path);

    /*---------------------------------------------------------
      Record 1: number of blocks
    ---------------------------------------------------------*/
    const auto block_count = static_cast<std::size_t>(reader.readScalar<std::int32_t>());
    
    /*---------------------------------------------------------
      Record 2: dimensions for each block (ni,nj,nk)
    ---------------------------------------------------------*/
    auto dimension_bytes = reader.readRecord();
    const std::size_t expected_bytes = block_count * 3 * sizeof(std::int32_t);
    if (dimension_bytes.size() != expected_bytes) {
        throw std::runtime_error("Unexpected dimension record size in Plot3D file.");
    }

    std::vector<std::array<int, 3>> dims(block_count);
    const auto* dims_ptr = reinterpret_cast<const std::int32_t*>(dimension_bytes.data());
    for (std::size_t block = 0; block < block_count; ++block) {
        dims[block][0] = static_cast<int>(dims_ptr[block * 3 + 0]);
        dims[block][1] = static_cast<int>(dims_ptr[block * 3 + 1]);
        dims[block][2] = static_cast<int>(dims_ptr[block * 3 + 2]);
    }

    /*---------------------------------------------------------
      Read coordinate records, one per block
    ---------------------------------------------------------*/
    std::vector<Plot3DZone> zones;
    zones.reserve(block_count);

    for (std::size_t block = 0; block < block_count; ++block) {
        const int ni = dims[block][0];
        const int nj = dims[block][1];
        const int nk = dims[block][2];
        const std::size_t point_count =
            static_cast<std::size_t>(ni) * nj * nk;
        const std::size_t value_count = point_count * 3;
        const std::size_t expected = value_count * sizeof(double);

        auto buffer = reader.readRecord();
        if (buffer.size() != expected) {
            throw std::runtime_error("Unexpected coordinate record size in Plot3D file.");
        }

        /*-----------------------------------------------------
          Unpack coordinates:
            [x1..xN, y1..yN, z1..zN]
        -----------------------------------------------------*/
        std::vector<double> raw(value_count);
        std::memcpy(raw.data(), buffer.data(), expected);

        Plot3DZone zone;
        zone.idx = static_cast<int>(block + 1);
        zone.name = "block-" + std::to_string(block + 1);
        zone.vtxSize = {static_cast<long long>(ni),
                        static_cast<long long>(nj),
                        static_cast<long long>(nk)};

        zone.x.resize(point_count);
        zone.y.resize(point_count);
        zone.z.resize(point_count);

        const double* x_src = raw.data();
        const double* y_src = raw.data() + point_count;
        const double* z_src = raw.data() + point_count * 2;
        
        for (std::size_t idx = 0; idx < point_count; ++idx) {
            zone.x[idx] = x_src[idx];
            zone.y[idx] = y_src[idx];
            zone.z[idx] = z_src[idx];
        }

        zones.push_back(std::move(zone));
    }

    return zones;
}

} // namespace bcdef
