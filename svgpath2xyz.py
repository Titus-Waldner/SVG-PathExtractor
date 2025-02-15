#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SVG Path Sampler

This module parses an SVG path definition (the d attribute of a <path> element)
and returns a dense list of (x, y, z) coordinates along the drawn path.
The third coordinate (z) is used as a pen flag: 0 means "pen down" (draw)
and 1 means "pen up" (lift, do not draw). The sampling density (i.e. the number
of sample points per segment) is configurable via the class constructor.

Based on SVGPATH2MPL code by Nezar Abdennur.
"""

from __future__ import division, print_function
from math import sin, cos, sqrt, degrees, radians, acos
import re
import numpy as np

# --- Regular expressions and constants ---
COMMANDS = set('MmZzLlHhVvCcSsQqTtAa')
UPPERCASE = set('MZLHVCSQTA')
COMMAND_RE = re.compile(r"([MmZzLlHhVvCcSsQqTtAa])")
FLOAT_RE = re.compile(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")

# -----------------------------------------------------------
# Functions for arc conversion (endpoint to center parameterization)
# -----------------------------------------------------------
def endpoint_to_center(start, radius, rotation, large, sweep, end):
    cosr = cos(radians(rotation))
    sinr = sin(radians(rotation))
    dx = (start.real - end.real) / 2
    dy = (start.imag - end.imag) / 2
    x1prim = cosr * dx + sinr * dy
    y1prim = -sinr * dx + cosr * dy
    x1prim_sq = x1prim * x1prim
    y1prim_sq = y1prim * y1prim

    rx = abs(radius.real)
    ry = abs(radius.imag)
    rx_sq = rx * rx
    ry_sq = ry * ry

    radius_scale = (x1prim_sq / rx_sq) + (y1prim_sq / ry_sq)
    if radius_scale > 1:
        radius_scale = sqrt(radius_scale)
        rx *= radius_scale
        ry *= radius_scale
        rx_sq = rx * rx
        ry_sq = ry * ry
    radius = rx + ry * 1j

    t1 = rx_sq * y1prim_sq
    t2 = ry_sq * x1prim_sq
    c = sqrt(abs((rx_sq * ry_sq - t1 - t2) / (t1 + t2)))
    if large == sweep:
        c = -c
    cxprim = c * rx * y1prim / ry
    cyprim = -c * ry * x1prim / rx

    center = complex(
        (cosr * cxprim - sinr * cyprim) + ((start.real + end.real) / 2),
        (sinr * cxprim + cosr * cyprim) + ((start.imag + end.imag) / 2)
    )

    ux = (x1prim - cxprim) / rx
    uy = (y1prim - cyprim) / ry
    vx = (-x1prim - cxprim) / rx
    vy = (-y1prim - cyprim) / ry
    n = sqrt(ux * ux + uy * uy)
    theta = degrees(acos(ux / n))
    if uy < 0:
        theta = -theta
    theta = theta % 360

    n = sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy))
    p = ux * vx + uy * vy
    d = p / n
    d = np.clip(d, -1.0, 1.0)
    delta = degrees(acos(d))
    if (ux * vy - uy * vx) < 0:
        delta = -delta
    delta = delta % 360
    if not sweep:
        if delta > 0:
            delta -= 360

    return radius, center, theta, theta + delta

def sample_arc(start, radius, rotation, large, sweep, end, density):
    r_complex, center, theta1, theta2 = endpoint_to_center(start, radius, rotation, large, sweep, end)
    theta1_rad = radians(theta1)
    theta2_rad = radians(theta2)
    theta_values = np.linspace(theta1_rad, theta2_rad, density)
    rx = abs(r_complex.real)
    ry = abs(r_complex.imag)
    rot_rad = radians(rotation)
    points = []
    for t in theta_values:
        x_prim = rx * cos(t)
        y_prim = ry * sin(t)
        x = cos(rot_rad) * x_prim - sin(rot_rad) * y_prim
        y = sin(rot_rad) * x_prim + cos(rot_rad) * y_prim
        points.append(center + complex(x, y))
    return points

# -----------------------------------------------------------
# Sampling functions for line and Bézier segments
# -----------------------------------------------------------
def sample_line(p0, p1, density):
    return [p0 + (p1 - p0) * t for t in np.linspace(0, 1, density)]

def sample_quadratic_bezier(p0, p1, p2, density):
    return [((1-t)**2) * p0 + 2 * (1-t) * t * p1 + (t**2) * p2 for t in np.linspace(0, 1, density)]

def sample_cubic_bezier(p0, p1, p2, p3, density):
    return [((1-t)**3) * p0 + 3 * ((1-t)**2) * t * p1 + 3 * (1-t) * (t**2) * p2 + (t**3) * p3
            for t in np.linspace(0, 1, density)]

# -----------------------------------------------------------
# Parsing functions (adapted from SVGPATH2MPL)
# -----------------------------------------------------------
def _tokenize_path(pathdef):
    for x in COMMAND_RE.split(pathdef):
        if x in COMMANDS:
            yield x
        for token in FLOAT_RE.findall(x):
            yield token

def _next_pos(elements):
    return float(elements.pop()) + float(elements.pop()) * 1j

def _parse_path(pathdef, current_pos):
    """Generator that yields (command, verts) tuples from the path definition string.
       Each verts is a list of (x, y) tuples.
    """
    elements = list(_tokenize_path(pathdef))
    elements.reverse()
    start_pos = None
    command = None

    while elements:
        if elements[-1] in COMMANDS:
            last_command = command
            command = elements.pop()
            absolute = command in UPPERCASE
            command = command.upper()
        else:
            if command is None:
                raise ValueError("Implicit command not allowed in path: " + pathdef)
            last_command = command

        if command == 'M':
            pos = _next_pos(elements)
            current_pos = pos if absolute else current_pos + pos
            start_pos = current_pos
            yield ('M', [(current_pos.real, current_pos.imag)])
            command = 'L'
        elif command == 'Z':
            if current_pos != start_pos:
                yield ('L', [(start_pos.real, start_pos.imag)])
            yield ('Z', [(start_pos.real, start_pos.imag)])
            current_pos = start_pos
            start_pos = None
            command = None
        elif command == 'L':
            pos = _next_pos(elements)
            pos = pos if absolute else current_pos + pos
            yield ('L', [(pos.real, pos.imag)])
            current_pos = pos
        elif command == 'H':
            x = float(elements.pop())
            pos = complex(x, current_pos.imag) if absolute else complex(current_pos.real + x, current_pos.imag)
            yield ('L', [(pos.real, pos.imag)])
            current_pos = pos
        elif command == 'V':
            y = float(elements.pop())
            pos = complex(current_pos.real, y) if absolute else complex(current_pos.real, current_pos.imag + y)
            yield ('L', [(pos.real, pos.imag)])
            current_pos = pos
        elif command == 'C':
            control1 = _next_pos(elements)
            control2 = _next_pos(elements)
            end = _next_pos(elements)
            if not absolute:
                control1 += current_pos
                control2 += current_pos
                end += current_pos
            verts = [(control1.real, control1.imag),
                     (control2.real, control2.imag),
                     (end.real, end.imag)]
            yield ('C', verts)
            current_pos = end
        elif command == 'S':
            if last_command not in ('C', 'S'):
                control1 = current_pos
            else:
                control1 = current_pos + (current_pos - prev_control2)
            control2 = _next_pos(elements)
            end = _next_pos(elements)
            if not absolute:
                control2 += current_pos
                end += current_pos
            verts = [(control1.real, control1.imag),
                     (control2.real, control2.imag),
                     (end.real, end.imag)]
            yield ('C', verts)
            prev_control2 = control2
            current_pos = end
        elif command == 'Q':
            control = _next_pos(elements)
            end = _next_pos(elements)
            if not absolute:
                control += current_pos
                end += current_pos
            verts = [(control.real, control.imag),
                     (end.real, end.imag)]
            yield ('Q', verts)
            current_pos = end
        elif command == 'T':
            if last_command not in ('Q', 'T'):
                control = current_pos
            else:
                control = current_pos + (current_pos - prev_control)
            end = _next_pos(elements)
            if not absolute:
                end += current_pos
            verts = [(control.real, control.imag),
                     (end.real, end.imag)]
            yield ('Q', verts)
            prev_control = control
            current_pos = end
        elif command == 'A':
            arc_radius = _next_pos(elements)
            rotation = float(elements.pop())
            large = bool(float(elements.pop()))
            sweep = bool(float(elements.pop()))
            end = _next_pos(elements)
            if not absolute:
                end += current_pos
            yield ('A', [(arc_radius.real, arc_radius.imag), rotation, large, sweep, (end.real, end.imag), (current_pos.real, current_pos.imag)])
            current_pos = end
        else:
            raise ValueError("Unknown command: " + command)

def sample_path(pathdef, current_pos=0+0j, density=20):
    """
    Parse the SVG path string and return a list of (x, y, z) tuples.
    z is 0 when the pen is down (draw) and 1 when the pen is up (lift).
    """
    segments = list(_parse_path(pathdef, current_pos))
    sampled = []  # list of tuples: (complex_point, pen_state)
    last_point = None
    subpath_start = None

    for seg in segments:
        command = seg[0]
        verts = seg[1]

        if command == 'M':
            pt = complex(verts[0][0], verts[0][1])
            # MOVETO: start a new subpath with pen lifted.
            sampled.append((pt, 1))
            last_point = pt
            subpath_start = pt
        elif command == 'L':
            pt = complex(verts[0][0], verts[0][1])
            if last_point is not None:
                pts = sample_line(last_point, pt, density)
                # All points on a drawn line: pen down (0)
                # (Omit the first point to avoid duplicate.)
                for p in pts[1:]:
                    sampled.append((p, 0))
            else:
                sampled.append((pt, 0))
            last_point = pt
        elif command == 'C':
            if last_point is None:
                continue
            p0 = last_point
            p1 = complex(verts[0][0], verts[0][1])
            p2 = complex(verts[1][0], verts[1][1])
            p3 = complex(verts[2][0], verts[2][1])
            pts = sample_cubic_bezier(p0, p1, p2, p3, density)
            for p in pts[1:]:
                sampled.append((p, 0))
            last_point = p3
        elif command == 'Q':
            if last_point is None:
                continue
            p0 = last_point
            p1 = complex(verts[0][0], verts[0][1])
            p2 = complex(verts[1][0], verts[1][1])
            pts = sample_quadratic_bezier(p0, p1, p2, density)
            for p in pts[1:]:
                sampled.append((p, 0))
            last_point = p2
        elif command == 'A':
            if last_point is None:
                continue
            arc_radius = complex(verts[0][0], verts[0][1])
            rotation = verts[1]
            large = verts[2]
            sweep = verts[3]
            end_pt = complex(verts[4][0], verts[4][1])
            pts = sample_arc(last_point, arc_radius, rotation, large, sweep, end_pt, density)
            for p in pts[1:]:
                sampled.append((p, 0))
            last_point = end_pt
        elif command == 'Z':
            # For closepath, we do not want to draw the connecting line.
            # Instead, we "lift" the pen at the end of the subpath.
            if last_point is not None and subpath_start is not None:
                # Optionally, you could sample a transition here.
                sampled.append((subpath_start, 1))
                last_point = subpath_start
            subpath_start = None
        else:
            # Ignore other commands.
            pass

    # Return as list of (x, y, z) tuples.
    return [(pt.real, pt.imag, pen) for pt, pen in sampled]

# -----------------------------------------------------------
# The SVGPathSampler class (wrapper around sample_path)
# -----------------------------------------------------------
class SVGPathSampler:
    def __init__(self, density=20):
        """
        density: number of sample points per segment.
        """
        self.density = density

    def sample_path(self, pathdef, current_pos=0+0j):
        return sample_path(pathdef, current_pos, self.density)

# -----------------------------------------------------------
# Example usage (for testing the module directly)
# -----------------------------------------------------------
if __name__ == '__main__':
    example_path = (
        "M 100 100 L 200 100 C 250 100 250 200 200 200 "
        "L 100 200 Z "
        "M 300 300 A 50 50 0 0 1 350 350"
    )
    sampler = SVGPathSampler(density=30)
    pts = sampler.sample_path(example_path)
    with open("sampled_points.txt", "w") as f:
        for x, y, z in pts:
            f.write("{:.2f} {:.2f} {:.0f}\n".format(x, y, z))
    print("Sampled", len(pts), "points along the SVG path.")
