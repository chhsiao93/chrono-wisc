// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2024 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Wei Hu, Huzaifa Mustafa Unjhawala
// =============================================================================
// Cratering validation problem involving spherical impactors with different
// densities falling from different heights (zero velocity) on CRM soil.
//
// Reference solution:
// https://www.sciencedirect.com/science/article/pii/S0045782521003534?ref=pdf_download&fr=RR-2&rr=8c4472d7d99222ff
//
// =============================================================================

#include <cassert>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/assets/ChVisualSystem.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"

#include "chrono_fsi/sph/ChFsiSystemSPH.h"
#include "chrono_fsi/sph/ChFsiFluidSystemSPH.h"

#ifdef CHRONO_VSG
    #include "chrono_fsi/sph/visualization/ChFsiVisualizationVSG.h"
#endif

#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/cxxopts/ChCLI.h"

#include <yaml-cpp/yaml.h>

using namespace chrono;
using namespace chrono::fsi;
using namespace chrono::fsi::sph;

// -----------------------------------------------------------------------------
// Note: sphere_radius is now a parameter that can be set via command-line
// -----------------------------------------------------------------------------

// Helper function to parse YAML array [x, y, z] to ChVector3d
ChVector3d YamlArrayToVector3d(const YAML::Node& node) {
    if (!node.IsSequence() || node.size() != 3) {
        std::cerr << "Error: Expected YAML array with 3 elements for ChVector3d" << std::endl;
        return ChVector3d(0, 0, 0);
    }
    return ChVector3d(node[0].as<double>(), node[1].as<double>(), node[2].as<double>());
}

#ifdef CHRONO_VSG
class MarkerPositionVisibilityCallback : public ChFsiVisualizationVSG::MarkerVisibilityCallback {
  public:
    MarkerPositionVisibilityCallback() {}
    virtual bool get(unsigned int n) const override { return pos[n].y > 0; }
};
#endif

// -----------------------------------------------------------------------------

// Function to handle CLI arguments
bool GetProblemSpecs(int argc,
                     char** argv,
                     double& t_end,
                     int& ps_freq,
                     double& output_fps,
                     double& sphere_density,
                     double& sphere_radius,
                     double& Hdrop,
                     double& initial_spacing,
                     double& d0_multiplier,
                     double& time_step,
                     ChVector3d& fluid_box_center,
                     ChVector3d& fluid_box_half,
                     bool& randomize_fluid_center,
                     ChVector3d& fluid_center_min,
                     ChVector3d& fluid_center_max,
                     bool& randomize_fluid_half,
                     ChVector3d& fluid_half_min,
                     ChVector3d& fluid_half_max,
                     bool& randomize_fluid_vel,
                     ChVector3d& fluid_vel,
                     ChVector3d& fluid_vel_min,
                     ChVector3d& fluid_vel_max,
                     bool& randomize_rigid_pos,
                     ChVector3d& rigid_pos,
                     ChVector3d& rigid_pos_min,
                     ChVector3d& rigid_pos_max,
                     bool& randomize_rigid_vel,
                     ChVector3d& rigid_vel,
                     ChVector3d& rigid_vel_min,
                     ChVector3d& rigid_vel_max,
                     int& random_seed,
                     std::string& output_folder,
                     std::string& boundary_type,
                     std::string& viscosity_type,
                     std::string& kernel_type) {
    ChCLI cli(argv[0], "FSI Cratering Demo - GNS Dataset Generation");

    // Primary configuration source: YAML file
    cli.AddOption<std::string>("Config", "config_file", "YAML configuration file path (required)", "");

    // Per-trajectory overrides (set by shell script for each trajectory)
    cli.AddOption<std::string>("Output", "output_folder", "Output folder name", "");
    cli.AddOption<double>("Geometry", "sphere_density", "Sphere density (overrides YAML)",
                          std::to_string(sphere_density));
    cli.AddOption<double>("Geometry", "sphere_radius", "Sphere radius (overrides YAML)", std::to_string(sphere_radius));
    cli.AddOption<int>("Random", "random_seed", "Random seed (0=use time)", std::to_string(random_seed));

    if (!cli.Parse(argc, argv))
        return false;

    // Load configuration from YAML file (required)
    std::string config_file = cli.GetAsType<std::string>("config_file");
    if (config_file.empty()) {
        std::cerr << "Error: --config_file is required!" << std::endl;
        std::cerr << "Usage: " << argv[0] << " --config_file <path_to_config.yaml>" << std::endl;
        return false;
    }

    try {
        YAML::Node config = YAML::LoadFile(config_file);

        // ===== Simulation Parameters =====
        if (config["simulation"]["t_end"]) {
            t_end = config["simulation"]["t_end"].as<double>();
        }
        if (config["simulation"]["time_step"]) {
            time_step = config["simulation"]["time_step"].as<double>();
        }
        if (config["simulation"]["ps_freq"]) {
            ps_freq = config["simulation"]["ps_freq"].as<int>();
        }
        if (config["simulation"]["output_fps"]) {
            output_fps = config["simulation"]["output_fps"].as<double>();
        }

        // ===== Geometry Parameters =====
        if (config["geometry"]["initial_spacing"]) {
            initial_spacing = config["geometry"]["initial_spacing"].as<double>();
        }
        if (config["geometry"]["Hdrop"]) {
            Hdrop = config["geometry"]["Hdrop"].as<double>();
        }

        // Fluid box vectors
        if (config["geometry"]["fluid_box"]["center"]) {
            fluid_box_center = YamlArrayToVector3d(config["geometry"]["fluid_box"]["center"]);
        }
        if (config["geometry"]["fluid_box"]["half_dim"]) {
            fluid_box_half = YamlArrayToVector3d(config["geometry"]["fluid_box"]["half_dim"]);
        }

        // ===== Physics Parameters =====
        if (config["physics"]["boundary_type"]) {
            boundary_type = config["physics"]["boundary_type"].as<std::string>();
        }
        if (config["physics"]["viscosity_type"]) {
            viscosity_type = config["physics"]["viscosity_type"].as<std::string>();
        }
        if (config["physics"]["kernel_type"]) {
            kernel_type = config["physics"]["kernel_type"].as<std::string>();
        }
        if (config["physics"]["d0_multiplier"]) {
            d0_multiplier = config["physics"]["d0_multiplier"].as<double>();
        }

        // ===== Fluid Box Center Randomization =====
        if (config["geometry"]["fluid_box_center"]["randomize"]) {
            randomize_fluid_center = config["geometry"]["fluid_box_center"]["randomize"].as<bool>();
        }
        if (config["geometry"]["fluid_box_center"]["min"]) {
            fluid_center_min = YamlArrayToVector3d(config["geometry"]["fluid_box_center"]["min"]);
        }
        if (config["geometry"]["fluid_box_center"]["max"]) {
            fluid_center_max = YamlArrayToVector3d(config["geometry"]["fluid_box_center"]["max"]);
        }

        // ===== Fluid Box Half Dim Randomization =====
        if (config["geometry"]["fluid_box_half_dim"]["randomize"]) {
            randomize_fluid_half = config["geometry"]["fluid_box_half_dim"]["randomize"].as<bool>();
        }
        if (config["geometry"]["fluid_box_half_dim"]["min"]) {
            fluid_half_min = YamlArrayToVector3d(config["geometry"]["fluid_box_half_dim"]["min"]);
        }
        if (config["geometry"]["fluid_box_half_dim"]["max"]) {
            fluid_half_max = YamlArrayToVector3d(config["geometry"]["fluid_box_half_dim"]["max"]);
        }

        // ===== Fluid Velocity Randomization =====
        if (config["fluid_velocity"]["randomize"]) {
            randomize_fluid_vel = config["fluid_velocity"]["randomize"].as<bool>();
        }
        if (config["fluid_velocity"]["min"]) {
            fluid_vel_min = YamlArrayToVector3d(config["fluid_velocity"]["min"]);
        }
        if (config["fluid_velocity"]["max"]) {
            fluid_vel_max = YamlArrayToVector3d(config["fluid_velocity"]["max"]);
        }

        // ===== Rigid Position =====
        if (config["rigid_position"]["randomize"]) {
            randomize_rigid_pos = config["rigid_position"]["randomize"].as<bool>();
        }
        if (config["rigid_position"]["min"]) {
            rigid_pos_min = YamlArrayToVector3d(config["rigid_position"]["min"]);
        }
        if (config["rigid_position"]["max"]) {
            rigid_pos_max = YamlArrayToVector3d(config["rigid_position"]["max"]);
        }

        // ===== Rigid Velocity =====
        if (config["rigid_velocity"]["randomize"]) {
            randomize_rigid_vel = config["rigid_velocity"]["randomize"].as<bool>();
        }
        if (config["rigid_velocity"]["min"]) {
            rigid_vel_min = YamlArrayToVector3d(config["rigid_velocity"]["min"]);
        }
        if (config["rigid_velocity"]["max"]) {
            rigid_vel_max = YamlArrayToVector3d(config["rigid_velocity"]["max"]);
        }

        std::cout << "Loaded configuration from: " << config_file << std::endl;
    } catch (const YAML::Exception& e) {
        std::cerr << "Error parsing YAML config file: " << e.what() << std::endl;
        return false;
    } catch (const std::exception& e) {
        std::cerr << "Error loading config file: " << e.what() << std::endl;
        return false;
    }

    // ===== CLI Overrides (for per-trajectory randomization) =====
    // These override YAML values if provided by shell script
    sphere_density = cli.GetAsType<double>("sphere_density");
    sphere_radius = cli.GetAsType<double>("sphere_radius");
    random_seed = cli.GetAsType<int>("random_seed");
    output_folder = cli.GetAsType<std::string>("output_folder");

    return true;
}
//------------------------------------------------------------------
// Function to generate a UV sphere mesh and save to VTK file
//------------------------------------------------------------------
void WriteSphereVTK(const std::string& filename, std::shared_ptr<ChBody> body, double radius, int resolution = 16) {
    // Generate a sphere mesh with UV coordinates
    ChTriangleMeshConnected mesh;
    std::vector<ChVector3d>& vertices = mesh.GetCoordsVertices();
    std::vector<ChVector3i>& indices = mesh.GetIndicesVertexes();

    // Create vertices using spherical coordinates
    for (int i = 0; i <= resolution; i++) {
        double phi = CH_PI * i / resolution;
        for (int j = 0; j <= resolution; j++) {
            double theta = 2 * CH_PI * j / resolution;
            double x = radius * sin(phi) * cos(theta);
            double y = radius * sin(phi) * sin(theta);
            double z = radius * cos(phi);
            vertices.push_back(ChVector3d(x, y, z));
        }
    }

    // Create triangular faces
    for (int i = 0; i < resolution; i++) {
        for (int j = 0; j < resolution; j++) {
            int p1 = i * (resolution + 1) + j;
            int p2 = p1 + 1;
            int p3 = (i + 1) * (resolution + 1) + j;
            int p4 = p3 + 1;

            // Add two triangles for each grid cell
            indices.push_back(ChVector3i(p1, p2, p3));
            indices.push_back(ChVector3i(p2, p4, p3));
        }
    }

    // Now write the mesh to VTK file, transforming it to the body's position
    std::ofstream outf(filename);
    outf << "# vtk DataFile Version 2.0" << std::endl;
    outf << "Sphere VTK from simulation" << std::endl;
    outf << "ASCII" << std::endl;
    outf << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write vertices transformed by body frame
    ChFrame<> frame = body->GetFrameRefToAbs();
    outf << "POINTS " << vertices.size() << " float" << std::endl;
    for (const auto& v : vertices) {
        auto w = frame.TransformPointLocalToParent(v);
        outf << w.x() << " " << w.y() << " " << w.z() << std::endl;
    }

    // Write triangular cells
    int nf = static_cast<int>(indices.size());
    outf << "CELLS " << nf << " " << 4 * nf << std::endl;
    for (const auto& f : indices) {
        outf << "3 " << f.x() << " " << f.y() << " " << f.z() << std::endl;
    }

    // Write cell types (5 = VTK_TRIANGLE)
    outf << "CELL_TYPES " << nf << std::endl;
    for (int i = 0; i < nf; i++) {
        outf << "5" << std::endl;
    }
    outf.close();
}

int main(int argc, char* argv[]) {
    // Default values
    double t_end = 2.0;
    bool verbose = true;
    bool output = true;
    bool write_marker_files = true;
    double output_fps = 400;
    bool snapshots = false;
    int ps_freq = 1;
    double sphere_density = 700;
    double sphere_radius = 0.0125;  // Default sphere radius (can be changed via CLI)
    double Hdrop = 0.2;
    bool render = false;
    double render_fps = 400;
    std::string boundary_type = "adami";
    std::string viscosity_type = "artificial_unilateral";
    double initial_spacing = 0.0025;
    double d0_multiplier = 1.3;
    double time_step = 5e-5;
    std::string kernel_type = "cubic";

    // Fluid box parameters (defaults based on original code)
    ChVector3d fluid_box_center(0.0, 0.0, 0.025);     // Center of fluid box (m)
    ChVector3d fluid_box_half(0.048, 0.048, 0.0225);  // Half dimensions (m)

    // Fluid box center randomization parameters
    bool randomize_fluid_center = false;              // If true, randomize fluid box center
    ChVector3d fluid_center_min(-0.01, -0.01, 0.015);  // Min fluid box center (m)
    ChVector3d fluid_center_max(0.01, 0.01, 0.025);    // Max fluid box center (m)

    // Fluid box half_dim randomization parameters
    bool randomize_fluid_half = false;                // If true, randomize fluid box half dimensions
    ChVector3d fluid_half_min(0.015, 0.015, 0.012);   // Min fluid box half dimensions (m)
    ChVector3d fluid_half_max(0.025, 0.025, 0.020);   // Max fluid box half dimensions (m)

    // Fluid velocity parameters
    bool randomize_fluid_vel = false;                // If true, randomize fluid velocity
    ChVector3d fluid_vel(0.0, 0.0, 0.0);             // Initial fluid velocity (m/s)
    ChVector3d fluid_vel_min(-0.1, -0.1, -0.1);      // Min fluid velocity (m/s)
    ChVector3d fluid_vel_max(0.1, 0.1, 0.1);         // Max fluid velocity (m/s)

    // Rigid body position parameters
    bool randomize_rigid_pos = false;              // If true, randomize position
    ChVector3d rigid_pos(0.0, 0.0, 0.0);           // Initial position (m) - will be calculated if not randomized
    ChVector3d rigid_pos_min(-0.05, -0.03, 0.15);  // Min position (m)
    ChVector3d rigid_pos_max(0.05, 0.03, 0.20);    // Max position (m)

    // Rigid body velocity parameters
    bool randomize_rigid_vel = false;            // If true, randomize velocity
    ChVector3d rigid_vel(0.0, 0.0, 0.0);         // Initial velocity (m/s)
    ChVector3d rigid_vel_min(-0.5, -0.5, -2.0);  // Min velocity (m/s)
    ChVector3d rigid_vel_max(0.5, 0.5, -0.5);    // Max velocity (m/s)

    int random_seed = 0;  // Random seed (0 = use system time)

    // Output folder name
    std::string output_folder = "";  // Empty = auto-generate from parameters

    // Parse command-line arguments
    if (!GetProblemSpecs(argc, argv, t_end, ps_freq, output_fps, sphere_density, sphere_radius, Hdrop, initial_spacing,
                         d0_multiplier, time_step, fluid_box_center, fluid_box_half, randomize_fluid_center,
                         fluid_center_min, fluid_center_max, randomize_fluid_half, fluid_half_min, fluid_half_max,
                         randomize_fluid_vel, fluid_vel, fluid_vel_min, fluid_vel_max, randomize_rigid_pos, rigid_pos,
                         rigid_pos_min, rigid_pos_max, randomize_rigid_vel, rigid_vel, rigid_vel_min, rigid_vel_max,
                         random_seed, output_folder, boundary_type, viscosity_type, kernel_type)) {
        return 1;
    }

    // Initialize random number generator if needed
    if (randomize_fluid_center || randomize_fluid_half || randomize_fluid_vel || randomize_rigid_pos || randomize_rigid_vel) {
        if (random_seed == 0) {
            random_seed = static_cast<int>(std::time(nullptr));
        }
        std::srand(random_seed);
        std::cout << "Using random seed: " << random_seed << std::endl;

        // Randomize fluid box center if enabled
        if (randomize_fluid_center) {
            fluid_box_center =
                fluid_center_min + ChVector3d((std::rand() / (double)RAND_MAX) * (fluid_center_max.x() - fluid_center_min.x()),
                                              (std::rand() / (double)RAND_MAX) * (fluid_center_max.y() - fluid_center_min.y()),
                                              (std::rand() / (double)RAND_MAX) * (fluid_center_max.z() - fluid_center_min.z()));
            std::cout << "Randomized fluid_box_center: (" << fluid_box_center.x() << ", " << fluid_box_center.y() << ", " << fluid_box_center.z()
                      << ") m" << std::endl;
        }

        // Randomize fluid box half dimensions if enabled
        if (randomize_fluid_half) {
            fluid_box_half =
                fluid_half_min + ChVector3d((std::rand() / (double)RAND_MAX) * (fluid_half_max.x() - fluid_half_min.x()),
                                            (std::rand() / (double)RAND_MAX) * (fluid_half_max.y() - fluid_half_min.y()),
                                            (std::rand() / (double)RAND_MAX) * (fluid_half_max.z() - fluid_half_min.z()));
            std::cout << "Randomized fluid_box_half: (" << fluid_box_half.x() << ", " << fluid_box_half.y() << ", " << fluid_box_half.z()
                      << ") m" << std::endl;
        }

        // Randomize fluid velocity if enabled
        if (randomize_fluid_vel) {
            fluid_vel =
                fluid_vel_min + ChVector3d((std::rand() / (double)RAND_MAX) * (fluid_vel_max.x() - fluid_vel_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_vel_max.y() - fluid_vel_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_vel_max.z() - fluid_vel_min.z()));
            std::cout << "Randomized fluid_vel: (" << fluid_vel.x() << ", " << fluid_vel.y() << ", " << fluid_vel.z()
                      << ") m/s" << std::endl;
        }

        // Randomize position if enabled
        if (randomize_rigid_pos) {
            rigid_pos =
                rigid_pos_min + ChVector3d((std::rand() / (double)RAND_MAX) * (rigid_pos_max.x() - rigid_pos_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_pos_max.y() - rigid_pos_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_pos_max.z() - rigid_pos_min.z()));
            std::cout << "Randomized position: (" << rigid_pos.x() << ", " << rigid_pos.y() << ", " << rigid_pos.z()
                      << ") m" << std::endl;
        }

        // Randomize velocity if enabled
        if (randomize_rigid_vel) {
            rigid_vel =
                rigid_vel_min + ChVector3d((std::rand() / (double)RAND_MAX) * (rigid_vel_max.x() - rigid_vel_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_vel_max.y() - rigid_vel_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_vel_max.z() - rigid_vel_min.z()));
            std::cout << "Randomized velocity: (" << rigid_vel.x() << ", " << rigid_vel.y() << ", " << rigid_vel.z()
                      << ") m/s" << std::endl;
        }
    }

    // ============================================================
    // Validate geometry and resample if overlapping
    // ============================================================
    const int MAX_RESAMPLE_ATTEMPTS = 100;
    int resample_count = 0;
    bool valid_configuration = false;

    // Define boundary dimensions (from the original code)
    double bxDim = 0.1;
    double byDim = 0.1;
    double bzDim = 0.1;

    // Calculate fluid box bounds
    ChVector3d fluid_min = fluid_box_center - fluid_box_half;
    ChVector3d fluid_max = fluid_box_center + fluid_box_half;

    while (!valid_configuration && resample_count < MAX_RESAMPLE_ATTEMPTS) {
        bool has_overlap = false;

        // Calculate default position if needed (before checks)
        if (!randomize_rigid_pos && rigid_pos.z() == 0.0) {
            rigid_pos = ChVector3d(0.0, 0.0, fluid_box_center.z() * 2 + sphere_radius + 0.5 * initial_spacing);
        }

        // Check 1: Sphere vs Fluid box overlap
        ChVector3d sphere_min = rigid_pos - ChVector3d(sphere_radius, sphere_radius, sphere_radius);
        ChVector3d sphere_max = rigid_pos + ChVector3d(sphere_radius, sphere_radius, sphere_radius);

        bool sphere_fluid_overlap =
            !(sphere_max.x() < fluid_min.x() || sphere_min.x() > fluid_max.x() || sphere_max.y() < fluid_min.y() ||
              sphere_min.y() > fluid_max.y() || sphere_max.z() < fluid_min.z() || sphere_min.z() > fluid_max.z());

        if (sphere_fluid_overlap) {
            std::cout << "Warning: Sphere overlaps with fluid box!" << std::endl;
            has_overlap = true;
        }

        // Check 2: Sphere vs Boundary walls
        ChVector3d boundary_min(-bxDim / 2, -byDim / 2, 0.0);
        ChVector3d boundary_max(bxDim / 2, byDim / 2, bzDim);  // Top is open

        double safety_margin = initial_spacing;  // Safety margin to prevent collision

        bool sphere_boundary_overlap =
            (sphere_min.x() < boundary_min.x() + safety_margin || sphere_max.x() > boundary_max.x() - safety_margin ||
             sphere_min.y() < boundary_min.y() + safety_margin || sphere_max.y() > boundary_max.y() - safety_margin ||
             sphere_min.z() < boundary_min.z() + safety_margin);

        if (sphere_boundary_overlap) {
            std::cout << "Warning: Sphere too close to boundary walls!" << std::endl;
            has_overlap = true;
        }

        // Check 3: Fluid box vs Boundary walls
        bool fluid_boundary_overlap =
            (fluid_min.x() < boundary_min.x() + safety_margin || fluid_max.x() > boundary_max.x() - safety_margin ||
             fluid_min.y() < boundary_min.y() + safety_margin || fluid_max.y() > boundary_max.y() - safety_margin ||
             fluid_min.z() < boundary_min.z() + safety_margin);

        if (fluid_boundary_overlap) {
            std::cout << "Warning: Fluid box too close to boundary walls!" << std::endl;
            has_overlap = true;
        }

        // If valid, exit loop
        if (!has_overlap) {
            valid_configuration = true;
            if (resample_count > 0) {
                std::cout << "Valid configuration found after " << resample_count << " resampling attempts"
                          << std::endl;
            }
        } else {
            // Resample if randomization is enabled
            if (randomize_fluid_center || randomize_fluid_half || randomize_fluid_vel || randomize_rigid_pos || randomize_rigid_vel) {
                resample_count++;
                std::cout << "Resampling attempt " << resample_count << "..." << std::endl;

                // Resample fluid box center if randomized
                if (randomize_fluid_center) {
                    fluid_box_center = fluid_center_min +
                                ChVector3d((std::rand() / (double)RAND_MAX) * (fluid_center_max.x() - fluid_center_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_center_max.y() - fluid_center_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_center_max.z() - fluid_center_min.z()));
                }

                // Resample fluid box half dimensions if randomized
                if (randomize_fluid_half) {
                    fluid_box_half = fluid_half_min +
                                ChVector3d((std::rand() / (double)RAND_MAX) * (fluid_half_max.x() - fluid_half_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_half_max.y() - fluid_half_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_half_max.z() - fluid_half_min.z()));
                }

                // Resample fluid velocity if randomized
                if (randomize_fluid_vel) {
                    fluid_vel = fluid_vel_min +
                                ChVector3d((std::rand() / (double)RAND_MAX) * (fluid_vel_max.x() - fluid_vel_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_vel_max.y() - fluid_vel_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (fluid_vel_max.z() - fluid_vel_min.z()));
                }

                // Resample position if randomized
                if (randomize_rigid_pos) {
                    rigid_pos = rigid_pos_min +
                                ChVector3d((std::rand() / (double)RAND_MAX) * (rigid_pos_max.x() - rigid_pos_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_pos_max.y() - rigid_pos_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_pos_max.z() - rigid_pos_min.z()));
                }

                // Resample velocity if randomized
                if (randomize_rigid_vel) {
                    rigid_vel = rigid_vel_min +
                                ChVector3d((std::rand() / (double)RAND_MAX) * (rigid_vel_max.x() - rigid_vel_min.x()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_vel_max.y() - rigid_vel_min.y()),
                                           (std::rand() / (double)RAND_MAX) * (rigid_vel_max.z() - rigid_vel_min.z()));
                }
            } else {
                // Cannot resample without randomization enabled
                std::cerr << "Error: Geometry overlap detected, but randomization is disabled!" << std::endl;
                std::cerr << "Please adjust parameters manually or enable randomization." << std::endl;
                return 1;
            }
        }
    }

    // Check if we exceeded maximum attempts
    if (!valid_configuration) {
        std::cerr << "Error: Could not find valid configuration after " << MAX_RESAMPLE_ATTEMPTS << " attempts!"
                  << std::endl;
        std::cerr << "Please adjust the parameter ranges to allow valid configurations." << std::endl;
        return 1;
    }

    std::cout << "Geometry validation passed!" << std::endl;
    std::cout << "Final rigid body position: (" << rigid_pos.x() << ", " << rigid_pos.y() << ", " << rigid_pos.z()
              << ")" << std::endl;
    std::cout << "Final rigid body velocity: (" << rigid_vel.x() << ", " << rigid_vel.y() << ", " << rigid_vel.z()
              << ")" << std::endl;

    // Create a physics system
    ChSystemSMC sysMBS;
    sysMBS.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Create a fluid system
    ChFsiFluidSystemSPH sysSPH;
    // Create an FSI system
    ChFsiSystemSPH sysFSI(sysMBS, sysSPH);

    sysFSI.SetStepSizeCFD(time_step);
    sysFSI.SetStepsizeMBD(time_step);

    ChFsiFluidSystemSPH::ElasticMaterialProperties mat_props;
    mat_props.density = 1510;
    mat_props.Young_modulus = 2e6;
    mat_props.Poisson_ratio = 0.3;
    mat_props.mu_I0 = 0.04;
    mat_props.mu_fric_s = 0.3;
    mat_props.mu_fric_2 = 0.48;
    mat_props.average_diam = 0.001;
    mat_props.cohesion_coeff = 0;  // default

    sysSPH.SetElasticSPH(mat_props);

    ChFsiFluidSystemSPH::SPHParameters sph_params;
    sph_params.sph_method = SPHMethod::WCSPH;
    sph_params.initial_spacing = initial_spacing;
    sph_params.d0_multiplier = d0_multiplier;
    sph_params.artificial_viscosity = 0.01;
    sph_params.shifting_method = ShiftingMethod::PPST_XSPH;
    sph_params.shifting_xsph_eps = 0.5;
    sph_params.shifting_ppst_pull = 1.0;
    sph_params.shifting_ppst_push = 3.0;
    sph_params.kernel_threshold = 0.8;
    sph_params.num_proximity_search_steps = ps_freq;
    if (kernel_type == "cubic") {
        sph_params.kernel_type = KernelType::CUBIC_SPLINE;
    } else if (kernel_type == "wendland") {
        sph_params.kernel_type = KernelType::WENDLAND;
    } else {
        std::cerr << "Invalid kernel type: " << kernel_type << std::endl;
        return 1;
    }

    // Set boundary type
    if (boundary_type == "holmes") {
        sph_params.boundary_type = BoundaryType::HOLMES;
    } else {
        sph_params.boundary_type = BoundaryType::ADAMI;
    }

    // Set viscosity type
    if (viscosity_type == "artificial_bilateral") {
        sph_params.viscosity_type = ViscosityType::ARTIFICIAL_BILATERAL;
    } else {
        sph_params.viscosity_type = ViscosityType::ARTIFICIAL_UNILATERAL;
    }

    sysSPH.SetSPHParameters(sph_params);

    double g = 9.81;
    sysSPH.SetGravitationalAcceleration(ChVector3d(0, 0, -g));
    sysMBS.SetGravitationalAcceleration(sysSPH.GetGravitationalAcceleration());

    sysFSI.SetVerbose(verbose);
    sysSPH.SetNumProximitySearchSteps(ps_freq);

    // Dimension of the space domain (use values already declared in validation section)
    double fzDim = fluid_box_center.z() * 2;  // Calculate fzDim from center if needed

    // Set the periodic boundary condition
    double init_spacing = initial_spacing;
    ChVector3d cMin(-bxDim / 2 - 3 * init_spacing, -byDim / 2 - 3 * init_spacing, -bzDim * 1.2);
    ChVector3d cMax(bxDim / 2 + 3 * init_spacing, byDim / 2 + 3 * init_spacing,
                    (bzDim + sphere_radius + init_spacing) * 1.2);
    sysSPH.SetComputationalBoundaries(cMin, cMax, PeriodicSide::NONE);

    // Create SPH particle locations using user-specified box parameters
    chrono::utils::ChGridSampler<> sampler(init_spacing);
    std::vector<ChVector3d> points = sampler.SampleBox(fluid_box_center, fluid_box_half);

    std::cout << "Fluid box center: (" << fluid_box_center.x() << ", " << fluid_box_center.y() << ", "
              << fluid_box_center.z() << ")" << std::endl;
    std::cout << "Fluid box half dimensions: (" << fluid_box_half.x() << ", " << fluid_box_half.y() << ", "
              << fluid_box_half.z() << ")" << std::endl;
    std::cout << "Generated " << points.size() << " SPH particles" << std::endl;

    // Add SPH particles to the fluid system
    double gz = std::abs(sysSPH.GetGravitationalAcceleration().z());
    // each particle's initial pressure and density based on depth
    for (const auto& p : points) {
        double pre_ini = 0;     // assign homogeneous initial pressure
        double rho_ini = 1800;  // assign homogeneous initial density
        sysSPH.AddSPHParticle(p, rho_ini, pre_ini, sysSPH.GetViscosity(), fluid_vel);
    }

    std::cout << "Fluid initial velocity: (" << fluid_vel.x() << ", " << fluid_vel.y() << ", " << fluid_vel.z() << ") m/s" << std::endl;

    // Create MBD and BCE particles for the solid domain
    auto cmaterial = chrono_types::make_shared<ChContactMaterialSMC>();
    cmaterial->SetYoungModulus(1e8);
    cmaterial->SetFriction(0.3f);
    cmaterial->SetRestitution(0.05f);
    cmaterial->SetAdhesion(0);

    // Create a container
    auto box = chrono_types::make_shared<ChBody>();
    box->SetPos(ChVector3d(0.0, 0.0, 0.0));
    box->SetRot(ChQuaternion<>(1, 0, 0, 0));
    box->SetFixed(true);
    // Add collision geometry for the container walls
    chrono::utils::AddBoxContainer(box, cmaterial,                                 //
                                   ChFrame<>(ChVector3d(0, 0, bzDim / 2), QUNIT),  //
                                   ChVector3d(bxDim, byDim, bzDim), 0.01,          //
                                   ChVector3i(2, 2, -1),                           //
                                   true);
    box->EnableCollision(true);
    sysMBS.AddBody(box);

    // Add BCE particles attached on the walls into FSI system
    sysSPH.AddBoxContainerBCE(box,                                            //
                              ChFrame<>(ChVector3d(0, 0, bzDim / 2), QUNIT),  //
                              ChVector3d(bxDim, byDim, bzDim),                //
                              ChVector3i(2, 2, -1));

    // Create a falling sphere
    double volume = ChSphere::GetVolume(sphere_radius);
    double mass = sphere_density * volume;
    auto inertia = mass * ChSphere::GetGyration(sphere_radius);

    // Calculate default velocity if not randomized/specified
    double impact_vel = 0.0;
    if (!randomize_rigid_vel && rigid_vel.z() == 0.0) {
        // Use default velocity based on Hdrop (impact velocity)
        impact_vel = std::sqrt(2 * Hdrop * g);
        rigid_vel = ChVector3d(0.0, 0.0, -impact_vel);  // Negative for downward
    }

    // Position and velocity are already set by validation loop above
    auto sphere = chrono_types::make_shared<ChBody>();
    sphere->SetPos(rigid_pos);
    sphere->SetPosDt(rigid_vel);
    sphere->SetMass(mass);
    sphere->SetInertia(inertia);

    chrono::utils::AddSphereGeometry(sphere.get(), cmaterial, sphere_radius, ChVector3d(0, 0, 0), QUNIT, true);
    sphere->EnableCollision(true);
    // sphere->GetCollisionModel()->SetSafeMargin(init_spacing);
    sysMBS.AddBody(sphere);

    sysFSI.AddFsiBody(sphere);
    sysSPH.AddSphereBCE(sphere, ChFrame<>(VNULL, QUNIT), sphere_radius, true, true);

    // Complete construction of the FSI system
    sysFSI.Initialize();

    // Output directories
    std::string out_dir;
    if (output || snapshots) {
        out_dir = GetChronoOutputPath() + "GNS_dataset/";
        if (!filesystem::create_directory(filesystem::path(out_dir))) {
            std::cerr << "Error creating directory " << out_dir << std::endl;
            return 1;
        }

        // Use user-specified folder name or auto-generate from parameters
        if (!output_folder.empty()) {
            out_dir = out_dir + output_folder;
        } else {
            std::stringstream ss;
            ss << viscosity_type << "_" << boundary_type;
            ss << "_ps" << ps_freq;
            ss << "_d" << sphere_density;
            ss << "_r" << sphere_radius;
            ss << "_h" << Hdrop;
            ss << "_s" << init_spacing;
            // Add velocity info to directory name if non-zero
            if (randomize_rigid_vel) {
                ss << "_vrand_seed" << random_seed;
            } else if (rigid_vel.x() != 0.0 || rigid_vel.y() != 0.0) {
                ss << "_vx" << rigid_vel.x() << "_vy" << rigid_vel.y();
            }
            out_dir = out_dir + ss.str();
        }
        if (!filesystem::create_directory(filesystem::path(out_dir))) {
            std::cerr << "Error creating directory " << out_dir << std::endl;
            return 1;
        }

        if (output) {
            if (!filesystem::create_directory(filesystem::path(out_dir + "/particles"))) {
                std::cerr << "Error creating directory " << out_dir + "/particles" << std::endl;
                return 1;
            }
            if (!filesystem::create_directory(filesystem::path(out_dir + "/fsi"))) {
                std::cerr << "Error creating directory " << out_dir + "/fsi" << std::endl;
                return 1;
            }
            if (!filesystem::create_directory(filesystem::path(out_dir + "/vtk"))) {
                std::cerr << "Error creating directory " << out_dir + "/vtk" << std::endl;
                return 1;
            }
        }

        if (snapshots) {
            if (!filesystem::create_directory(filesystem::path(out_dir + "/snapshots"))) {
                std::cerr << "Error creating directory " << out_dir + "/snapshots" << std::endl;
                return 1;
            }
        }
    }

    // Create a run-time visualizer
    std::shared_ptr<ChVisualSystem> vis;

#ifdef CHRONO_VSG
    if (render) {
        // FSI plugin
        auto col_callback = chrono_types::make_shared<ParticleVelocityColorCallback>(0, impact_vel / 2);

        auto visFSI = chrono_types::make_shared<ChFsiVisualizationVSG>(&sysFSI);
        visFSI->EnableFluidMarkers(true);
        visFSI->EnableBoundaryMarkers(true);
        visFSI->EnableRigidBodyMarkers(false);
        visFSI->SetSPHColorCallback(col_callback);
        visFSI->SetSPHVisibilityCallback(chrono_types::make_shared<MarkerPositionVisibilityCallback>());
        visFSI->SetBCEVisibilityCallback(chrono_types::make_shared<MarkerPositionVisibilityCallback>());

        // VSG visual system (attach visFSI as plugin)
        auto visVSG = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
        visVSG->AttachPlugin(visFSI);
        visVSG->AttachSystem(&sysMBS);
        visVSG->SetWindowTitle("Cratering");
        visVSG->SetWindowSize(1280, 800);
        visVSG->SetWindowPosition(100, 100);
        visVSG->AddCamera(ChVector3d(0, -3 * byDim, 0.75 * bzDim), ChVector3d(0, 0, 0.75 * bzDim));
        visVSG->SetLightIntensity(0.9f);
        visVSG->SetLightDirection(-CH_PI_2, CH_PI / 6);

        visVSG->Initialize();
        vis = visVSG;
    }
#else
    render = false;
#endif

    // Start the simulation
    double time = 0.0;
    int sim_frame = 0;
    int out_frame = 0;
    int render_frame = 0;
    double dT = sysFSI.GetStepSizeCFD();

    double rtf_average = 0.0;
    unsigned int rtf_count = 0;

    std::string out_file = out_dir + "/sphere_penetration_depth.txt";
    std::ofstream ofile(out_file, std::ios::trunc);

    ChTimer timer;
    timer.start();
    while (time < t_end) {
        if (output && time >= out_frame / output_fps) {
            if (write_marker_files) {
                sysSPH.SaveParticleData(out_dir + "/particles");
                sysSPH.SaveSolidData(out_dir + "/fsi", time);
                std::stringstream vtk_filename;
                vtk_filename << out_dir << "/vtk/sphere_" << std::setw(5) << std::setfill('0') << out_frame << ".vtk";
                WriteSphereVTK(vtk_filename.str(), sphere, sphere_radius);
                // std::cout << " -- Wrote " << vtk_filename.str() << std::endl;
            }
            // std::cout << " -- Output frame " << out_frame << " at t = " << time << std::endl;
            out_frame++;
        }

        if (render && time >= render_frame / render_fps) {
            if (!vis->Run())
                break;
            vis->Render();

            if (snapshots) {
                std::cout << " -- Snapshot frame " << render_frame << " at t = " << time << std::endl;
                std::ostringstream filename;
                filename << out_dir << "/snapshots/" << std::setw(5) << std::setfill('0') << render_frame << ".jpg";
                vis->WriteImageToFile(filename.str());
            }

            render_frame++;
        }

        // Write penetration depth to file
        double d_pen = fzDim + sphere_radius + 0.5 * init_spacing - sphere->GetPos().z();
        ofile << time << " " << d_pen << " " << sphere->GetPos().x() << " " << sphere->GetPos().y() << " "
              << sphere->GetPos().z() << " " << sphere->GetPosDt().x() << " " << sphere->GetPosDt().y() << " "
              << sphere->GetPosDt().z() << std::endl;
        // Advance simulation for one timestep for all systems
        sysFSI.DoStepDynamics(dT);
        double rtf = sysFSI.GetRtf();
        rtf_average = (rtf_average * rtf_count + rtf) / (rtf_count + 1);
        rtf_count++;
        time += dT;

        sim_frame++;
    }
    timer.stop();
    std::cout << "End Time: " << t_end << std::endl;
    std::cout << "\nSimulation time: " << timer() << " seconds\n" << std::endl;
    std::cout << "Average RTF: " << rtf_average << std::endl;
    std::cout << "Simulation finished" << std::endl;

    // Write runtime information and all parameters to a file
    if (output) {
        std::ofstream runtime_file(out_dir + "/runtime.txt");

        // Runtime information
        runtime_file << "=== Runtime Information ===" << std::endl;
        runtime_file << "Runtime: " << timer() << " seconds" << std::endl;
        runtime_file << "Simulation time: " << time << " seconds" << std::endl;
        runtime_file << "Average RTF: " << rtf_average << std::endl;
        runtime_file << std::endl;

        // Simulation parameters
        runtime_file << "=== Simulation Parameters ===" << std::endl;
        runtime_file << "t_end: " << t_end << std::endl;
        runtime_file << "time_step: " << time_step << std::endl;
        runtime_file << "ps_freq: " << ps_freq << std::endl;
        runtime_file << std::endl;

        // Geometry parameters
        runtime_file << "=== Geometry Parameters ===" << std::endl;
        runtime_file << "sphere_density: " << sphere_density << " kg/m^3" << std::endl;
        runtime_file << "sphere_radius: " << sphere_radius << " m" << std::endl;
        runtime_file << "Hdrop: " << Hdrop << " m" << std::endl;
        runtime_file << "initial_spacing: " << initial_spacing << " m" << std::endl;
        runtime_file << std::endl;

        // Fluid box parameters
        runtime_file << "=== Fluid Box Parameters ===" << std::endl;
        runtime_file << "randomize_fluid_center: " << (randomize_fluid_center ? "true" : "false") << std::endl;
        runtime_file << "fluid_box_center: (" << fluid_box_center.x() << ", " << fluid_box_center.y() << ", "
                     << fluid_box_center.z() << ") m" << std::endl;
        if (randomize_fluid_center) {
            runtime_file << "fluid_center_min: (" << fluid_center_min.x() << ", " << fluid_center_min.y() << ", "
                         << fluid_center_min.z() << ") m" << std::endl;
            runtime_file << "fluid_center_max: (" << fluid_center_max.x() << ", " << fluid_center_max.y() << ", "
                         << fluid_center_max.z() << ") m" << std::endl;
        }
        runtime_file << "randomize_fluid_half: " << (randomize_fluid_half ? "true" : "false") << std::endl;
        runtime_file << "fluid_box_half: (" << fluid_box_half.x() << ", " << fluid_box_half.y() << ", "
                     << fluid_box_half.z() << ") m" << std::endl;
        if (randomize_fluid_half) {
            runtime_file << "fluid_half_min: (" << fluid_half_min.x() << ", " << fluid_half_min.y() << ", "
                         << fluid_half_min.z() << ") m" << std::endl;
            runtime_file << "fluid_half_max: (" << fluid_half_max.x() << ", " << fluid_half_max.y() << ", "
                         << fluid_half_max.z() << ") m" << std::endl;
        }
        runtime_file << "randomize_fluid_vel: " << (randomize_fluid_vel ? "true" : "false") << std::endl;
        runtime_file << "fluid_vel: (" << fluid_vel.x() << ", " << fluid_vel.y() << ", " << fluid_vel.z() << ") m/s"
                     << std::endl;
        if (randomize_fluid_vel) {
            runtime_file << "fluid_vel_min: (" << fluid_vel_min.x() << ", " << fluid_vel_min.y() << ", "
                         << fluid_vel_min.z() << ") m/s" << std::endl;
            runtime_file << "fluid_vel_max: (" << fluid_vel_max.x() << ", " << fluid_vel_max.y() << ", "
                         << fluid_vel_max.z() << ") m/s" << std::endl;
        }
        runtime_file << "num_sph_particles: " << points.size() << std::endl;
        runtime_file << std::endl;

        // Rigid body position parameters
        runtime_file << "=== Rigid Body Position ===" << std::endl;
        runtime_file << "randomize_rigid_pos: " << (randomize_rigid_pos ? "true" : "false") << std::endl;
        runtime_file << "rigid_pos: (" << rigid_pos.x() << ", " << rigid_pos.y() << ", " << rigid_pos.z() << ") m"
                     << std::endl;
        if (randomize_rigid_pos) {
            runtime_file << "rigid_pos_min: (" << rigid_pos_min.x() << ", " << rigid_pos_min.y() << ", "
                         << rigid_pos_min.z() << ") m" << std::endl;
            runtime_file << "rigid_pos_max: (" << rigid_pos_max.x() << ", " << rigid_pos_max.y() << ", "
                         << rigid_pos_max.z() << ") m" << std::endl;
        }
        runtime_file << std::endl;

        // Rigid body velocity parameters
        runtime_file << "=== Rigid Body Velocity ===" << std::endl;
        runtime_file << "randomize_rigid_vel: " << (randomize_rigid_vel ? "true" : "false") << std::endl;
        runtime_file << "rigid_vel: (" << rigid_vel.x() << ", " << rigid_vel.y() << ", " << rigid_vel.z() << ") m/s"
                     << std::endl;
        if (randomize_rigid_vel) {
            runtime_file << "rigid_vel_min: (" << rigid_vel_min.x() << ", " << rigid_vel_min.y() << ", "
                         << rigid_vel_min.z() << ") m/s" << std::endl;
            runtime_file << "rigid_vel_max: (" << rigid_vel_max.x() << ", " << rigid_vel_max.y() << ", "
                         << rigid_vel_max.z() << ") m/s" << std::endl;
        }
        if (randomize_rigid_pos || randomize_rigid_vel) {
            runtime_file << "random_seed: " << random_seed << std::endl;
        }
        runtime_file << std::endl;

        // Physics parameters
        runtime_file << "=== Physics Parameters ===" << std::endl;
        runtime_file << "boundary_type: " << boundary_type << std::endl;
        runtime_file << "viscosity_type: " << viscosity_type << std::endl;
        runtime_file << "kernel_type: " << kernel_type << std::endl;
        runtime_file << "d0_multiplier: " << d0_multiplier << std::endl;
        runtime_file << std::endl;

        runtime_file << "Simulation finished" << std::endl;
        runtime_file.close();
    }

    return 0;
}