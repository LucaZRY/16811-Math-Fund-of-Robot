#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt


def findConvexHull(points):
    """
    Andrew's monotone chain convex hull.
    Input: iterable of 2D points (N x 2).
    Output: hull vertices in CCW order (M x 2).
    """
    pts = np.asarray(points, dtype=float)
    if pts.shape[0] < 3:
        return pts.copy()

    # sort lexicographically (x, then y)
    pts = pts[np.lexsort((pts[:, 1], pts[:, 0]))]

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # build lower hull
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(tuple(p))

    # build upper hull
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(tuple(p))

    # concatenate lower and upper, removing last point of each (duplicate)
    hull = np.array(lower[:-1] + upper[:-1], dtype=float)
    return hull


def plot_polygon(points, hull):
    points = np.asarray(points, float)
    hull = np.asarray(hull, float)

    plt.figure()
    plt.scatter(points[:, 0], points[:, 1], color="blue", label="points")

    if hull.shape[0] >= 1:
        xs = np.append(hull[:, 0], hull[0, 0])
        ys = np.append(hull[:, 1], hull[0, 1])
        plt.plot(xs, ys, "r-", linewidth=2, label="convex hull")

    plt.axis("equal")
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.title("Convex Hull (Part 1)")
    plt.tight_layout()
    plt.show()


def load_points_from_file(fname):
    pts = []
    with open(fname, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # allow "x y" or "x, y"
            for sep in [",", " "]:
                line = line.replace(sep, " ")
            toks = [t for t in line.split() if t]
            if len(toks) != 2:
                raise ValueError(f"Bad line in {fname}: {line}")
            x, y = map(float, toks)
            pts.append([x, y])
    return np.array(pts, float)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        points = load_points_from_file(sys.argv[1])
    else:
        # random example if no file provided
        np.random.seed(0)
        points = np.random.rand(20, 2) * 3.0

    hull = findConvexHull(points)
    print("Hull vertices (in order):")
    print(hull)
    plot_polygon(points, hull)
