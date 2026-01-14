/*─────────────────────────────────────────────────────────────
  File: src/bc_define_main.cpp
  CLI entry point and high-level flow
─────────────────────────────────────────────────────────────*/
#include "bc_define_internal.hpp"
#include "logger.hpp"
#include "mesh_utils.hpp"

#include <filesystem>
#include <iostream>

namespace {

void usage()
{
    std::cerr << "Usage: bc_define mesh.cgns bcdef.input "
                 "[--overwrite] [--autowall] [--autofarfield] [--threads N]\n"
                 "       bc_define mesh.cgns [bcdef.input] --bench-io [--bench-iter N]\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        usage();
        return 1;
    }

    const std::string meshPath = argv[1];
    std::string bcdefPath;

    bool overwrite = false;
    bool autowall = false;
    bool autofarfield = false;
    bool bench_io = false;
    int threads = 0;
    int bench_iters = 3;

    int argi = 2;
    if (argi < argc && argv[argi][0] != '-') {
        bcdefPath = argv[argi];
        ++argi;
    }

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

    const bool is_plot3d = fs::bc::is_plot3d_mesh(meshPath);

    std::filesystem::path logPath = std::filesystem::path(meshPath).parent_path() /
                                    "bc_define.log";
    fs::Logger log(logPath.string());
    log.info("BC_define started");
    log.info("Mesh     : " + meshPath);
    if (!bcdefPath.empty())
        log.info("BC file  : " + bcdefPath);

    if (autowall && autofarfield) {
        std::cerr << "Error: --autowall and --autofarfield are mutually exclusive\n";
        return 1;
    }

    if (bench_io && is_plot3d) {
        std::cerr << "Error: --bench-io is only supported for CGNS meshes\n";
        return 1;
    }

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
        if (is_plot3d) {
            if (bcdefPath.empty()) {
                usage();
                return 1;
            }

            auto zones = fs::read_plot3d(meshPath);
            log.info("Opened Plot3D mesh with " + std::to_string(zones.size()) + " zone(s)");

            auto result = fs::ConnectivityDetector::run_plot3d(
                zones, log, fs::GEOM_TOL, threads);
            auto bc_specs = fs::bc::parse_bcdef(bcdefPath);

            auto bcs = fs::bc::build_plot3d_bcs(
                zones, result.boundary, bc_specs, autowall, autofarfield);
            std::filesystem::path cwd = std::filesystem::current_path();
            fs::bc::write_plot3d_1to1s((cwd / "1to1s").string(), zones, result.connections);
            fs::bc::write_plot3d_bcs((cwd / "bcs").string(), zones, bcs);

            log.info("Finished OK.  Exiting.");
            return 0;
        }

        fs::Mesh mesh;
        mesh.open(meshPath, /*modify=*/!bench_io);

        if (bench_io) {
            fs::bc::run_io_benchmark(mesh, log, bench_iters);
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

        if (overwrite)
            fs::purge_BC_and_connectivity(mesh, log);

        std::vector<fs::BoundaryPatch> boundary_patches =
            fs::ConnectivityDetector::run(mesh, log, fs::GEOM_TOL, overwrite, threads);

        if (!boundary_patches.empty()) {
            auto bc_specs = fs::bc::parse_bcdef(bcdefPath);
            fs::bc::write_boundary_conditions(mesh, boundary_patches, bc_specs,
                                              autowall, autofarfield);
        }

        if (overwrite)
            fs::prune_unused_families(mesh, log);

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
