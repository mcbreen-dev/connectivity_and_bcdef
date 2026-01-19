/*─────────────────────────────────────────────────────────────
  File: src/bc_define_main.cpp

  Command-line entry point for boundary-condition definition
  and multizone connectivity processing.

  This file:
  - Parses CLI arguments and validates combinations
  - Dispatches to either CGNS or Plot3D workflows
  - Orchestrates connectivity detection, BC definition, and writing
  - Handles optional destructive overwrite and I/O benchmarking

  All geometry processing, connectivity detection, and BC writing
  are delegated to other modules; this file contains no mesh math.
─────────────────────────────────────────────────────────────*/
#include "bc_define_internal.hpp"
#include "logger.hpp"
#include "cgns_cleanup.hpp"

#include <filesystem>
#include <iostream>

namespace {

/*-------------------------------------------------------------
  Print CLI usage and exit conditions.

  The CLI supports two mutually exclusive operating modes:
  1) Normal BC + connectivity definition
  2) CGNS coordinate I/O benchmarking (--bench-io)
-------------------------------------------------------------*/
void usage()
{
    std::cerr << "Usage: bc_define mesh.cgns bcdef.input "
                 "[--overwrite] [--autowall] [--autofarfield] [--threads N]\n"
                 "       bc_define mesh.cgns [bcdef.input] --bench-io [--bench-iter N]\n";
}

} // anonymous namespace


/*=====================================================================
  Program entry point.

  Control flow overview:

  1) Parse CLI arguments
  2) Determine mesh format (CGNS vs Plot3D)
  3) Initialize logging
  4) Dispatch to:
     - Plot3D connectivity + BC file generation, OR
     - CGNS connectivity detection and in-file BC writing
=====================================================================*/
int main(int argc, char** argv)
{
    /*---------------------------------------------------------
      Basic argument validation: at least a mesh path is required
    ---------------------------------------------------------*/
    if (argc < 2) {
        usage();
        return 1;
    }

    const std::string meshPath = argv[1];
    std::string bcdefPath;

    /*---------------------------------------------------------
      CLI flags with defaults.

      overwrite      : delete existing BC/connectivity nodes first
      autowall       : assign BCWall to unspecified boundary faces
      autofarfield   : assign BCFarfield to unspecified boundary faces
      bench_io       : run coordinate I/O benchmark instead of BC logic
      threads        : worker thread count (0 → auto)
      bench_iters    : benchmark iterations per read
    ---------------------------------------------------------*/
    bool overwrite = false;
    bool autowall = false;
    bool autofarfield = false;
    bool bench_io = false;
    int threads = 0;
    int bench_iters = 3;

    /*---------------------------------------------------------
      Optional positional bcdef file:
      - Must appear immediately after mesh path
      - Must not begin with '-'
    ---------------------------------------------------------*/
    int argi = 2;
    if (argi < argc && argv[argi][0] != '-') {
        bcdefPath = argv[argi];
        ++argi;
    }


    /*---------------------------------------------------------
      Parse remaining flags
    ---------------------------------------------------------*/
    for (int i = argi; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--overwrite") overwrite = true;
        else if (flag == "--autowall") autowall = true;
        else if (flag == "--autofarfield") autofarfield = true;
        else if (flag == "--bench-io") bench_io = true;
        else if (flag == "--bench-iter") {
            if (i + 1 >= argc) {
                usage();
                return 1;
            }
            bench_iters = std::stoi(argv[++i]);
        }
        else if (flag == "--threads") {
            if (i + 1 >= argc) {
                usage();
                return 1;
            }
            threads = std::stoi(argv[++i]);
        } else if (flag.rfind("--threads=", 0) == 0) {
            threads = std::stoi(flag.substr(std::string("--threads=").size()));
        }
        else {
            usage();
            return 1;
        }
    }

    /*---------------------------------------------------------
      Mesh format detection:
      Plot3D meshes are identified purely by filename extension.
    ---------------------------------------------------------*/
    const bool is_plot3d = bcdef::boundary::is_plot3d_mesh(meshPath);

    /*---------------------------------------------------------
      Logger initialization.
      The log file is written next to the mesh file.
    ---------------------------------------------------------*/
    std::filesystem::path logPath = std::filesystem::path(meshPath).parent_path() /
                                    "bc_define.log";
    bcdef::Logger log(logPath.string());
    log.info("BC_define started");
    log.info("Mesh     : " + meshPath);
    if (!bcdefPath.empty())
        log.info("BC file  : " + bcdefPath);

    /*---------------------------------------------------------
      Flag validation
    ---------------------------------------------------------*/
    if (autowall && autofarfield) {
        std::cerr << "Error: --autowall and --autofarfield are mutually exclusive\n";
        return 1;
    }

    if (bench_io && is_plot3d) {
        std::cerr << "Error: --bench-io is only supported for CGNS meshes\n";
        return 1;
    }

    /*---------------------------------------------------------
      Log effective operating mode
    ---------------------------------------------------------*/
    if (bench_io) {
        log.info("Flags    : bench_io=yes iters=" + std::to_string(bench_iters));
    } else {
        std::string flags = "Flags    : overwrite=" + std::string(overwrite ? "yes" : "no") +
                            " autowall=" + std::string(autowall ? "yes" : "no") +
                            " autofarfield=" + std::string(autofarfield ? "yes" : "no") +
                            " threads=" + std::to_string(threads);
        log.info(flags);
    }

    try {
        /*=====================================================
          Plot3D workflow:
          - Read mesh
          - Detect connectivity
          - Write text-based 1to1s and BC files
        =====================================================*/
        if (is_plot3d) {
            if (bcdefPath.empty()) {
                usage();
                return 1;
            }

            auto zones = bcdef::read_plot3d(meshPath);
            log.info("Opened Plot3D mesh with " + std::to_string(zones.size()) + " zone(s)");

            auto result = bcdef::ConnectivityDetector::run_plot3d(
                zones, log, bcdef::GEOM_TOL, threads);
            auto bc_specs = bcdef::boundary::parse_bcdef(bcdefPath);

            auto bcs = bcdef::boundary::build_plot3d_bcs(
                zones, result.boundary, bc_specs, autowall, autofarfield);
            std::filesystem::path cwd = std::filesystem::current_path();
            bcdef::boundary::write_plot3d_1to1s((cwd / "1to1s").string(), zones, result.connections);
            bcdef::boundary::write_plot3d_bcs((cwd / "bcs").string(), zones, bcs);

            log.info("Finished OK.  Exiting.");
            return 0;
        }

        /*=====================================================
          CGNS workflow
        =====================================================*/
        bcdef::Mesh mesh;
        mesh.open(meshPath, /*modify=*/!bench_io);

        /*-----------------------------------------------------
          Optional I/O benchmark mode:
          No BCs or connectivity are modified.
        -----------------------------------------------------*/
        if (bench_io) {
            bcdef::boundary::run_io_benchmark(mesh, log, bench_iters);
            mesh.close();
            log.info("Finished OK.  Exiting.");
            return 0;
        }

        if (bcdefPath.empty()) {
            usage();
            return 1;
        }

        if (!overwrite) {
            log.info("Opened mesh with " + std::to_string(mesh.zones().size()) +
                     " zone(s)");
        }

        /*-----------------------------------------------------
          Optional destructive cleanup:
          Removes existing BC_t and GridConnectivity_t nodes.
        -----------------------------------------------------*/
        if (overwrite)
            bcdef::purge_BC_and_connectivity(mesh, log);

/*
CGNS-SIDS design intent:
A given zone boundary segment should be covered by either (a) a BoundaryCondition_t entry, or
(b) a multizone interface connectivity entry (e.g., GridConnectivity1to1_t), but not both.
This program treats detected 1-to-1 matches as connectivity and assigns BCs only to remaining
unmatched boundary patches; --overwrite removes existing BC/connectivity nodes before rewriting.

CGNS User Guide. CGNS/SIDS - Standard Interface Data Structures/8:Multizone Interface Connectivity
https://cgns.org/standard/SIDS/multizone.html
*/

        std::vector<bcdef::BoundaryPatch> boundary_patches =
            bcdef::ConnectivityDetector::run(mesh, log, bcdef::GEOM_TOL, overwrite, threads);

        /*-----------------------------------------------------
          Write BCs only if uncovered boundary patches exist
        -----------------------------------------------------*/
        if (!boundary_patches.empty()) {
            auto bc_specs = bcdef::boundary::parse_bcdef(bcdefPath);
            bcdef::boundary::write_boundary_conditions(mesh, boundary_patches, bc_specs,
                                              autowall, autofarfield);
        }

        /*-----------------------------------------------------
          Remove unused CGNS Family_t nodes after overwrite
        -----------------------------------------------------*/
        if (overwrite)
            bcdef::prune_unused_families(mesh, log);

        mesh.close();
        log.info("Finished OK.  Exiting.");
        return 0;
    }
    catch (const std::exception& ex) {
        log.error(ex.what());
        std::cout << "Fatal: " << ex.what() << '\n';
        return 2;
    }
}
