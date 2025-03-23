"""
svg_mapper.py

This script uses your SVGParser (from the module 'command') to generate an output file from an SVG.
It then loads SVG path strings from that output file and uses an instance of SVGPathSampler (from
the module 'svgpath2xy') to sample the drawn path into a dense list of (x, y, z) points.
Points with z = 0 are drawn (pen down) and points with z = 1 indicate a pen lift (do not draw).
The script then plots the drawing by drawing only the segments where the pen is down.
"""


import re
import matplotlib.pyplot as plt
import numpy as np
from command import SVGParser       
from svgpath2xyz import SVGPathSampler

def load_svg_path_strings(command_file):
    """
    Load SVG path definition strings from a file.
    Lines starting with 'Path:' are expected to have the path data enclosed in quotes.
    Returns a list of SVG path strings.
    """
    path_strings = []
    # This regex expects lines like:  Path: "M 100 100 L 200 100 ... z"
    path_pattern = re.compile(r'\s*\s*"(.*)"')
    with open(command_file, 'r') as f:
        for line in f:
            m = path_pattern.search(line)
            if m:
                d = m.group(1)
                path_strings.append(d)
    return path_strings

def main(input = "input.svg", output = "output.txt"):
    # 1. Use your SVGParser to generate the output file from the SVG.
    parser = SVGParser(input, output)
    parser.parse_svg()
    
    # 2. Load SVG path definition strings from output.txt.
    svg_paths = load_svg_path_strings("output.txt")
    if not svg_paths:
        print("No SVG paths found in output.txt")
        return
    
    # 3. Create an instance of SVGPathSampler with desired density.
    sampler = SVGPathSampler(density=100)
    
    # 4. Create a matplotlib figure.
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # loop over each SVG path string.
    for d in svg_paths:
        try:
            pts = sampler.sample_path(d)
        except Exception as e:
            print("Error sampling path:", d, "\n", e)
            continue
        if not pts:
            continue
        
        # write the sampled points to a file.
        with open("sampled_points.txt", "w") as f:
            for x, y, z in pts:
                f.write("{:.2f} {:.2f} {:.0f}\n".format(x, y, z))
            print("Sampled", len(pts), "points along the SVG path.")
        
        # process the points into segments.
        segments = []
        current_segment = []
        pen_up_positions = []  # To store pen-up positions
        
        for (x, y, pen) in pts:
            if pen == 1:
                # pen up
                if current_segment:
                    segments.append(current_segment)
                    current_segment = []
                pen_up_positions.append((x, y))
            else:
                # Pen down
                current_segment.append((x, y))
        # Add any remaining segment.
        if current_segment:
            segments.append(current_segment)
        
        #  assign each segment a color from the viridis colormap.
        total_segments = len(segments)
        for idx, seg in enumerate(segments):
            xs, ys = zip(*seg)
            # normalize the index to get a color gradient.
            color_factor = idx / max(total_segments - 1, 1)
            color = plt.cm.viridis(color_factor)
            ax.plot(xs, ys, linestyle='-', color=color, linewidth=2, label=f"Segment {idx+1}" if total_segments > 1 else "Path")
        
        # pen-up positions as red circle markers.
        if pen_up_positions:
            pen_x, pen_y = zip(*pen_up_positions)
            ax.scatter(pen_x, pen_y, color='red', marker='o', s=50, zorder=5, label="Pen Lift")
    
    # plot 4 visualizaitons
    ax.set_xlim(0, 800)
    ax.set_ylim(0, 800)
    ax.set_aspect('equal')
    ax.set_title("Visualization of Sampled SVG Paths")
    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")
    ax.grid(True)
    plt.gca().invert_yaxis()  # invert y-axis to match SVG coordinate system
    
    # legend
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys())
    
    plt.show()

if __name__ == "__main__":
    main("smile_face.svg")
