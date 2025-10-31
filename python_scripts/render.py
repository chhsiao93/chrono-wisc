import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import imageio.v2 as imageio
import open3d as o3d
from matplotlib import cm
import argparse
import sys

# ============================================================
# Parse command-line arguments
# ============================================================
def parse_args():
    parser = argparse.ArgumentParser(
        description='Render FSI simulation data to video',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required argument
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to simulation data directory (containing particles/ folder)')

    # Optional rendering settings
    parser.add_argument('--width', type=int, default=1280,
                        help='Render width in pixels')
    parser.add_argument('--height', type=int, default=720,
                        help='Render height in pixels')
    parser.add_argument('--point_size', type=float, default=4.0,
                        help='Point size for particles')
    parser.add_argument('--fps', type=int, default=10,
                        help='Frames per second for output video')
    parser.add_argument('--format', type=str, default='mp4', choices=['mp4', 'gif'],
                        help='Output format: mp4 or gif')
    parser.add_argument('--jpeg_quality', type=int, default=85,
                        help='JPEG quality (0-100) for frame compression')
    parser.add_argument('--output_name', type=str, default='simulation',
                        help='Output video filename (without extension)')

    return parser.parse_args()

# Parse arguments
args = parse_args()

# Set parameters from arguments
output_dir = args.data_path
RENDER_WIDTH = args.width
RENDER_HEIGHT = args.height
POINT_SIZE = args.point_size
FPS = args.fps
OUTPUT_FORMAT = args.format
JPEG_QUALITY = args.jpeg_quality
OUTPUT_NAME = args.output_name

# Validate data path
if not os.path.exists(output_dir):
    print(f"Error: Data path does not exist: {output_dir}")
    sys.exit(1)

particles_dir = os.path.join(output_dir, "particles")
if not os.path.exists(particles_dir):
    print(f"Error: Particles directory not found: {particles_dir}")
    sys.exit(1)

print(f"Loading data from: {output_dir}")
print(f"Rendering settings: {RENDER_WIDTH}x{RENDER_HEIGHT}, {FPS} fps, format={OUTPUT_FORMAT}")

# Read boundary box for setting plot limits
boundary_file = f"{output_dir}/particles/boundary0.csv"
boundary_data = pd.read_csv(boundary_file)
boundary_pos = boundary_data[['x', 'y', 'z']].to_numpy()
min_bounds = np.min(boundary_pos, axis=0)
max_bounds = np.max(boundary_pos, axis=0)

# Find all fluid CSV files and sort them numerically
def extract_frame_number(filename):
    """Extract frame number from filename like 'fluid123.csv' or 'rigidBCE123.csv'"""
    import re
    match = re.search(r'(\d+)\.csv$', filename)
    return int(match.group(1)) if match else 0

fluid_files = sorted(glob.glob(f"{output_dir}/particles/fluid*.csv"), key=extract_frame_number)
rigid_files = sorted(glob.glob(f"{output_dir}/particles/rigidBCE*.csv"), key=extract_frame_number)
# print(f"Found {len(fluid_files)} fluid timesteps and {len(rigid_files)} rigid timesteps")

# Create render directory inside the data path
render_dir = os.path.join(output_dir, "render")
os.makedirs(render_dir, exist_ok=True)
print(f"Output directory: {render_dir}")

# Store image filenames for GIF creation
image_files = []

# Loop through all frames
num_frames = min(len(fluid_files), len(rigid_files))

# Find |U| min and max for color scaling
u_min, u_max = float('inf'), float('-inf')
for frame in range(num_frames):
    try:
        fluid_data = pd.read_csv(fluid_files[frame])
        u_values = fluid_data['|U|'].to_numpy()
        u_min = min(u_min, np.min(u_values))
        u_max = max(u_max, np.max(u_values))
    except Exception as e:
        print(f"Error reading {fluid_files[frame]}: {e}")
        continue

cmap = plt.get_cmap('viridis')
norm = plt.Normalize(vmin=u_min, vmax=u_max)

print(f"Rendering {num_frames} frames...")

# Setup Open3D offscreen renderer (works on headless systems)
renderer = o3d.visualization.rendering.OffscreenRenderer(RENDER_WIDTH, RENDER_HEIGHT)

# Setup rendering properties
renderer.scene.set_background([1.0, 1.0, 1.0, 1.0])  # White background

# Setup material for fluid point clouds
fluid_mat = o3d.visualization.rendering.MaterialRecord()
fluid_mat.shader = 'defaultUnlit'
fluid_mat.point_size = POINT_SIZE
# Setup material for rigid point clouds
rigid_mat = o3d.visualization.rendering.MaterialRecord()
rigid_mat.shader = 'defaultUnlit'
rigid_mat.point_size = POINT_SIZE * 2.0  # Rigid particles slightly larger


# Calculate camera parameters based on scene bounds
center = (min_bounds + max_bounds) / 2
scene_size = np.linalg.norm(max_bounds - min_bounds)
camera_distance = scene_size * 0.6

for frame in range(num_frames):
    try:
        # Read frame data
        fluid_data = pd.read_csv(fluid_files[frame])
        rigid_data = pd.read_csv(rigid_files[frame])

        fluid_pos = fluid_data[['x', 'y', 'z', '|U|']].to_numpy()
        rigid_pos = rigid_data[['x', 'y', 'z']].to_numpy()

        if frame == 0:
            print(f"Fluid particles: {fluid_pos.shape[0]}, Rigid BCE particles: {rigid_pos.shape[0]}")

        # Create point clouds
        fluid_pcd = o3d.geometry.PointCloud()
        fluid_pcd.points = o3d.utility.Vector3dVector(fluid_pos[:, :3])

        # Color fluid particles by velocity magnitude
        velocities = fluid_pos[:, 3]
        colors = cm.viridis(norm(velocities))[:, :3]  # Get RGB, drop alpha
        fluid_pcd.colors = o3d.utility.Vector3dVector(colors)

        # Create rigid particle point cloud (black color)
        rigid_pcd = o3d.geometry.PointCloud()
        rigid_pcd.points = o3d.utility.Vector3dVector(rigid_pos)
        rigid_pcd.paint_uniform_color([0, 0, 0])  # Black

        # Clear scene and add geometries
        renderer.scene.clear_geometry()
        renderer.scene.add_geometry("fluid", fluid_pcd, fluid_mat)
        renderer.scene.add_geometry("rigid", rigid_pcd, rigid_mat)

        # Setup camera on first frame
        if frame == 0:
            renderer.setup_camera(60.0, # field of view
                                center, # lookat
                                center + [camera_distance, camera_distance, camera_distance], # camera position
                                [0, 0, 1]) # up vector

        # Render to image
        img = renderer.render_to_image()

        # Save frame as JPEG (much smaller than PNG)
        frame_path = os.path.join(render_dir, f"frame_{frame:04d}.jpg")
        o3d.io.write_image(frame_path, img, quality=JPEG_QUALITY)
        image_files.append(frame_path)
    except Exception as e:
        print(f"Error rendering frame {frame}: {e}")
        continue

    if (frame + 1) % 10 == 0:
        print(f"  Rendered {frame + 1}/{num_frames} frames")

if OUTPUT_FORMAT == 'mp4':
    print("Creating MP4 video...")
    # Create MP4 using imageio-ffmpeg (much smaller than GIF)
    output_path = os.path.join(render_dir, f"{OUTPUT_NAME}.mp4")
    writer = imageio.get_writer(output_path, fps=FPS, codec='libx264',
                                 quality=8, pixelformat='yuv420p')

    for filename in image_files:
        writer.append_data(imageio.imread(filename))

    writer.close()
    print(f"MP4 video saved to {output_path}")

elif OUTPUT_FORMAT == 'gif':
    print("Creating GIF...")
    # Read all images and create GIF
    images = []
    for filename in image_files:
        images.append(imageio.imread(filename))

    # Save as optimized GIF
    gif_path = os.path.join(render_dir, f"{OUTPUT_NAME}.gif")
    duration = 1000 / FPS  # Convert FPS to milliseconds per frame
    imageio.mimsave(gif_path, images, duration=duration, loop=0,
                    subrectangles=True)
    print(f"GIF saved to {gif_path}")

print(f"Total frames: {num_frames}")
print(f"Rendering complete!")
