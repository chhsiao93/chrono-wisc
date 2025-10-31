#!/bin/bash

# ============================================================
# GNS Dataset Generation Script (Config-Based Version)
# ============================================================
# This script reads parameters from a YAML configuration file
# and generates multiple trajectories for GNS training
# ============================================================

# Configuration file path (default: config.yaml)
CONFIG_FILE="${1:-config.yaml}"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not found: $CONFIG_FILE"
    echo "Usage: $0 [config_file.yaml]"
    echo ""
    echo "Please create a config file based on config_template.yaml"
    exit 1
fi

echo "============================================================"
echo "GNS Dataset Generation - Config-Based"
echo "============================================================"
echo "Loading configuration from: $CONFIG_FILE"
echo ""

# ============================================================
# Function to parse YAML file (simple key-value parser)
# ============================================================
parse_yaml() {
    local yaml_file=$1
    local prefix=$2
    local s='[[:space:]]*'
    local w='[a-zA-Z0-9_]*'
    local fs=$(echo @|tr @ '\034')

    sed -ne "s|^\($s\):|\1|" \
         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p" "$yaml_file" |
    awk -F$fs '{
        indent = length($1)/2;
        vname[indent] = $2;
        for (i in vname) {if (i > indent) {delete vname[i]}}
        if (length($3) > 0) {
            vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
            printf("%s%s=\"%s\"\n", "'$prefix'", vn $2, $3);
        }
    }'
}

# Load configuration
eval $(parse_yaml "$CONFIG_FILE" "cfg_")

# ============================================================
# Helper function to get config value with default
# ============================================================
get_config() {
    local var_name="cfg_$1"
    local default_value="$2"
    echo "${!var_name:-$default_value}"
}

# ============================================================
# Read only the minimal configuration needed for shell script
# (C++ reads everything else directly from the YAML file!)
# ============================================================

# Dataset generation settings
NUM_TRAJECTORIES=$(get_config "dataset_num_trajectories" "10")
OUTPUT_BASE_DIR=$(get_config "dataset_output_base_dir" "GNS_dataset")
RANDOM_SEED_BASE=$(get_config "dataset_random_seed_base" "0")

# Sphere randomization ranges (for sampling random values each trajectory)
SPHERE_DENSITY_MIN=$(get_config "geometry_sphere_density_min" "700")
SPHERE_DENSITY_MAX=$(get_config "geometry_sphere_density_max" "2200")
SPHERE_RADIUS_MIN=$(get_config "geometry_sphere_radius_min" "0.008")
SPHERE_RADIUS_MAX=$(get_config "geometry_sphere_radius_max" "0.015")

# Rendering parameters (for auto-render after simulation)
AUTO_RENDER=$(get_config "rendering_auto_render" "true")
RENDER_WIDTH=$(get_config "rendering_width" "1280")
RENDER_HEIGHT=$(get_config "rendering_height" "720")
RENDER_FPS=$(get_config "rendering_fps" "10")
RENDER_FORMAT=$(get_config "rendering_format" "mp4")
CLEANUP_FRAMES=$(get_config "rendering_cleanup_frames" "true")
PYTHON_VENV=$(get_config "rendering_python_venv" "")

RENDER_SCRIPT="../../python_scripts/render.py"

# ============================================================
# Print configuration summary
# ============================================================
echo "Configuration Summary:"
echo "  Config file: $CONFIG_FILE"
echo "  Trajectories: $NUM_TRAJECTORIES"
echo "  Output dir: DEMO_OUTPUT/$OUTPUT_BASE_DIR"
echo ""
echo "Randomization Ranges (shell script samples these):"
echo "  Sphere density: [$SPHERE_DENSITY_MIN, $SPHERE_DENSITY_MAX] kg/m³"
echo "  Sphere radius: [$SPHERE_RADIUS_MIN, $SPHERE_RADIUS_MAX] m"
echo ""
echo "Other Parameters (read by C++ from YAML):"
echo "  - Fluid box geometry"
echo "  - Rigid position ranges"
echo "  - Rigid velocity ranges"
echo "  - Physics settings"
echo ""
echo "Rendering: $AUTO_RENDER (${RENDER_WIDTH}x${RENDER_HEIGHT}, ${RENDER_FPS}fps, $RENDER_FORMAT)"
echo "============================================================"
echo ""

# Activate Python virtual environment if specified
if [ "$AUTO_RENDER" = "true" ] && [ -n "$PYTHON_VENV" ]; then
    if [ -f "$PYTHON_VENV/bin/activate" ]; then
        echo "Activating Python virtual environment: $PYTHON_VENV"
        source "$PYTHON_VENV/bin/activate"
        echo "Using Python: $(which python3)"
        echo ""
    else
        echo "Warning: Virtual environment not found at $PYTHON_VENV/bin/activate"
        echo "Continuing with system Python..."
        echo ""
    fi
fi

# ============================================================
# Random number generation functions
# ============================================================
random_float() {
    local min=$1
    local max=$2
    awk -v min="$min" -v max="$max" 'BEGIN{srand(); printf "%.6f", min+rand()*(max-min)}'
}

# ============================================================
# Generate trajectories
# ============================================================
for i in $(seq 1 $NUM_TRAJECTORIES); do
    # Generate random parameters
    sphere_radius=$(random_float $SPHERE_RADIUS_MIN $SPHERE_RADIUS_MAX)
    sphere_density=$(random_float $SPHERE_DENSITY_MIN $SPHERE_DENSITY_MAX)

    # Generate unique seed for this trajectory
    if [ "$RANDOM_SEED_BASE" -eq 0 ]; then
        traj_seed=$RANDOM
    else
        traj_seed=$((RANDOM_SEED_BASE + i))
    fi

    # Format trajectory name
    traj_name=$(printf "traj%03d" $i)

    echo "------------------------------------------------------------"
    echo "Generating trajectory $i/$NUM_TRAJECTORIES: $traj_name"
    echo "  sphere_radius:  $sphere_radius m"
    echo "  sphere_density: $sphere_density kg/m³"
    echo "  random_seed:    $traj_seed"
    echo "------------------------------------------------------------"

    # Build command - now much simpler with config file!
    cmd="./GNS_benchmark_Acc_FSI_Cratering"
    cmd="$cmd --config_file $CONFIG_FILE"
    cmd="$cmd --output_folder $traj_name"
    cmd="$cmd --sphere_density $sphere_density"
    cmd="$cmd --sphere_radius $sphere_radius"
    cmd="$cmd --random_seed $traj_seed"

    # Run simulation
    eval $cmd

    # Check if simulation was successful
    if [ $? -eq 0 ]; then
        echo "Simulation completed successfully!"

        # Auto-render if enabled
        if [ "$AUTO_RENDER" = "true" ]; then
            echo ""
            echo "Rendering trajectory $traj_name..."

            output_path="DEMO_OUTPUT/$OUTPUT_BASE_DIR/$traj_name"

            if [ -f "$RENDER_SCRIPT" ]; then
                python3 $RENDER_SCRIPT \
                    --data_path "$output_path" \
                    --width $RENDER_WIDTH \
                    --height $RENDER_HEIGHT \
                    --fps $RENDER_FPS \
                    --format $RENDER_FORMAT \
                    --output_name "simulation"

                if [ $? -eq 0 ]; then
                    echo "Rendering completed: $output_path/render/simulation.$RENDER_FORMAT"

                    # Clean up frame images if enabled
                    if [ "$CLEANUP_FRAMES" = "true" ]; then
                        echo "Cleaning up frame images..."
                        frame_count=$(ls -1 $output_path/render/frame_*.jpg 2>/dev/null | wc -l)
                        if [ $frame_count -gt 0 ]; then
                            rm -f $output_path/render/frame_*.jpg
                            echo "Deleted $frame_count frame images"
                        fi
                    fi
                else
                    echo "Warning: Rendering failed for $traj_name"
                fi
            else
                echo "Warning: Render script not found at $RENDER_SCRIPT"
            fi
        fi
    else
        echo "Error: Simulation failed for $traj_name"
    fi

    echo ""
done

echo "============================================================"
echo "All $NUM_TRAJECTORIES trajectories completed!"
echo "============================================================"
