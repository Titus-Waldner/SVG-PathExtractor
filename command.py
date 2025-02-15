import xml.etree.ElementTree as ET
import math

class SVGParser:
    def __init__(self, svg_file, output_file="command.txt"):
        self.svg_file = svg_file
        self.output_file = output_file
        self.namespace = {'svg': 'http://www.w3.org/2000/svg'}
    
    @staticmethod
    def parse_style(style_str):
        style = {}
        for item in style_str.split(';'):
            if ':' in item:
                key, value = item.split(':', 1)
                style[key.strip()] = value.strip()
        return style

    @staticmethod
    def parse_float(value, default=0.0):
        try:
            return float(value)
        except (ValueError, TypeError):
            return default

    @staticmethod
    def approximate_circle(cx, cy, r, segments=36):
        points = [(cx + r * math.cos(2 * math.pi * i / segments),
                   cy + r * math.sin(2 * math.pi * i / segments))
                  for i in range(segments)]
        points.append(points[0])
        return points

    @staticmethod
    def approximate_ellipse(cx, cy, rx, ry, segments=36):
        points = [(cx + rx * math.cos(2 * math.pi * i / segments),
                   cy + ry * math.sin(2 * math.pi * i / segments))
                  for i in range(segments)]
        points.append(points[0])
        return points

    @staticmethod
    def rect_to_lines(x, y, width, height):
        return [
            ((x, y), (x + width, y)),
            ((x + width, y), (x + width, y + height)),
            ((x + width, y + height), (x, y + height)),
            ((x, y + height), (x, y))
        ]

    @staticmethod
    def write_line_segments(f, points, shape_name="Line"):
        for (x1, y1), (x2, y2) in zip(points, points[1:]):
            f.write(f"{shape_name}: ({x1:.2f}, {y1:.2f}) -> ({x2:.2f}, {y2:.2f})\n")

    def parse_svg(self):
        tree = ET.parse(self.svg_file)
        root = tree.getroot()
        with open(self.output_file, "w") as f:
            
            for rect in root.findall(".//svg:rect", self.namespace):
                x = self.parse_float(rect.get("x", "0"))
                y = self.parse_float(rect.get("y", "0"))
                width = self.parse_float(rect.get("width", "0"))
                height = self.parse_float(rect.get("height", "0"))
                segments = self.rect_to_lines(x, y, width, height)
                self.write_line_segments(f, [seg[0] for seg in segments] + [segments[0][0]])

            for circle in root.findall(".//svg:circle", self.namespace):
                cx = self.parse_float(circle.get("cx", "0"))
                cy = self.parse_float(circle.get("cy", "0"))
                r = self.parse_float(circle.get("r", "0"))
                points = self.approximate_circle(cx, cy, r)
                self.write_line_segments(f, points)

            for ellipse in root.findall(".//svg:ellipse", self.namespace):
                cx = self.parse_float(ellipse.get("cx", "0"))
                cy = self.parse_float(ellipse.get("cy", "0"))
                rx = self.parse_float(ellipse.get("rx", "0"))
                ry = self.parse_float(ellipse.get("ry", "0"))
                points = self.approximate_ellipse(cx, cy, rx, ry)
                self.write_line_segments(f, points)

            for line in root.findall(".//svg:line", self.namespace):
                x1, y1, x2, y2 = line.get("x1"), line.get("y1"), line.get("x2"), line.get("y2")
                f.write(f"Line: ({x1}, {y1}) -> ({x2}, {y2})\n")

            for path in root.findall(".//svg:path", self.namespace):
                d = path.get("d", "")
                f.write(f'"{d}"\n')
        
        print(f"SVG commands extracted to {self.output_file}")

# Example usage:
# parser = SVGParser("input.svg")
# parser.parse_svg()
