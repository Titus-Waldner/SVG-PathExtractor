#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
from command import SVGParser       # Your existing SVGParser class
from svgpath2xyz import SVGPathSampler  # Our sampler module that returns (x,y,z) points

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

def main():
    # 1. Use your SVGParser to generate the output file from the SVG.
    parser = SVGParser("input.svg", "output.txt")
    parser.parse_svg()
    
    # 2. Load SVG path definition strings from output.txt.
    svg_paths = load_svg_path_strings("output.txt")
    if not svg_paths:
        print("No SVG paths found in output.txt")
        return
    
    # 3. Create an instance of SVGPathSampler with desired density.
    sampler = SVGPathSampler(density=100)
    
    # 4. For each SVG path string, sample to (x,y,z) points.
    # We'll split the points into segments based on the pen state.
    fig, ax = plt.subplots(figsize=(10, 10))
    
    for d in svg_paths:
        try:
            pts = sampler.sample_path(d)
        except Exception as e:
            print("Error sampling path:", d, "\n", e)
            continue
        if not pts:
            continue
        
        # Separate the sampled points into segments that should be drawn.
        segment_xs = []
        segment_ys = []
        # Iterate through points; when a point with pen=1 is encountered,
        # finish the current segment (if any) and then start a new segment.
        for (x, y, pen) in pts:
            if pen == 1:
                # Pen up: finish current segment if it exists.
                if segment_xs and segment_ys:
                    ax.plot(segment_xs, segment_ys, linestyle='-', color='blue')
                    segment_xs = []
                    segment_ys = []
                # Optionally, you could mark the pen-up position with a marker.
                # For now, we simply do not draw a line to it.
            else:
                # Pen down: add point to current segment.
                segment_xs.append(x)
                segment_ys.append(y)
        # Draw any remaining segment.
        if segment_xs and segment_ys:
            ax.plot(segment_xs, segment_ys, linestyle='-', color='blue')
    
    ax.set_xlim(0, 800)
    ax.set_ylim(0, 800)
    ax.set_aspect('equal')
    ax.set_title("Sampled (x,y,z) Points from SVG Paths (z=0 drawn, z=1 lifted)")
    plt.gca().invert_yaxis()  # Match the SVG coordinate system
    plt.show()

if __name__ == "__main__":
    main()
