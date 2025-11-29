#!/usr/bin/env python3
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from heapq import heappush, heappop

from q1 import findConvexHull


# ===========================
# Geometry utilities
# ===========================

def seg_intersect(a, b, c, d):
    """
    Proper segment intersection test.
    Returns True if segments ab and cd intersect in their interiors.
    Touching at endpoints is allowed, NOT considered an intersection.
    """

    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])

    o1 = orient(a, b, c)
    o2 = orient(a, b, d)
    o3 = orient(c, d, a)
    o4 = orient(c, d, b)

    # general case: strict crossing
    if o1 * o2 < 0 and o3 * o4 < 0:
        return True

    return False  # touching is allowed


def point_in_poly(point, poly):
    """
    Check if point strictly inside a convex polygon.

    Assumes polygon vertices are in CCW order.
    For a convex CCW polygon, a point is inside/on boundary if all
    cross products (edge -> point) are >= 0. Here we treat "strictly inside"
    for our collision check by using this as a heuristic on the midpoint.
    """
    x0, y0 = point
    inside = True
    for i in range(len(poly)):
        a = poly[i]
        b = poly[(i + 1) % len(poly)]
        cross = (b[0] - a[0]) * (y0 - a[1]) - (b[1] - a[1]) * (x0 - a[0])
        if cross < 0:  # outside for CCW polygon
            inside = False
            break
    return inside


def segment_clear(a, b, polygons):
    """
    Check if segment a–b is collision-free w.r.t. given convex polygons.

    Rules:
    - Segment may touch polygon boundaries (vertices/edges).
    - Segment may NOT cross the interior of any polygon.
    """
    for P in polygons:
        n = len(P)

        # 1) segment intersection with polygon edges (interior crossing)
        for i in range(n):
            c = P[i]
            d = P[(i + 1) % n]
            if seg_intersect(a, b, c, d):
                return False

        # 2) midpoint check for passing through polygon interior
        mid = (a + b) / 2.0
        if point_in_poly(mid, P):
            return False

    return True


# ===========================
# Visibility Graph
# ===========================

class VisibilityGraph:
    def __init__(self, polygons, start, goal):
        """
        polygons: list of convex polygons (each is (M_i x 2) numpy array)
        start, goal: 2D numpy arrays or lists
        """
        self.polygons = [np.array(p, float) for p in polygons]
        self.start = np.array(start, float)
        self.goal = np.array(goal, float)

        # Step 1: validity check (start/goal inside obstacle → no path)
        self.check_valid()

        # Step 2: build nodes = [start, goal] + all polygon vertices
        self.nodes = [self.start, self.goal]
        for P in self.polygons:
            for v in P:
                self.nodes.append(v)
        self.nodes = np.array(self.nodes)

        # Step 3: build edges (visibility graph)
        self.edges = {}
        self.build_edges()

        # Step 4: run Dijkstra
        self.path = self.dijkstra()

    def check_valid(self):
        """
        If start or goal lies strictly inside any polygon, we declare no path.
        """
        for P in self.polygons:
            if point_in_poly(self.start, P) or point_in_poly(self.goal, P):
                raise ValueError("Start or goal lies inside a polygon → no path.")

    def build_edges(self):
        n = len(self.nodes)
        for i in range(n):
            for j in range(i + 1, n):
                a, b = self.nodes[i], self.nodes[j]
                if segment_clear(a, b, self.polygons):
                    self.edges.setdefault(i, []).append(j)
                    self.edges.setdefault(j, []).append(i)

    def dijkstra(self):
        """
        Run Dijkstra on the visibility graph from node 0 (start) to node 1 (goal).
        Returns a list of 2D points (numpy arrays) along the shortest path,
        or None if no path exists.
        """
        start = 0
        goal = 1
        INF = float("inf")

        dist = {i: INF for i in range(len(self.nodes))}
        prev = {i: None for i in range(len(self.nodes))}
        dist[start] = 0.0

        pq = [(0.0, start)]
        while pq:
            d, u = heappop(pq)
            if d > dist[u]:
                continue
            if u == goal:
                break

            for v in self.edges.get(u, []):
                w = np.linalg.norm(self.nodes[u] - self.nodes[v])
                nd = dist[u] + w
                if nd < dist[v]:
                    dist[v] = nd
                    prev[v] = u
                    heappush(pq, (nd, v))

        if dist[goal] == INF:
            return None  # no path

        # reconstruct path indices
        path_indices = []
        cur = goal
        while cur is not None:
            path_indices.append(cur)
            cur = prev[cur]
        path_indices.reverse()

        return [self.nodes[i] for i in path_indices]

    def plot(self, title=None, save_path=None):
        """
        Plot environment, visibility edges, and shortest path.
        If save_path is provided, save the figure to that file and close it.
        Otherwise, show it interactively.
        """
        plt.figure()

        # polygons
        for P in self.polygons:
            xs = np.append(P[:, 0], P[0, 0])
            ys = np.append(P[:, 1], P[0, 1])
            plt.plot(xs, ys, "k-")

        # visibility edges
        for i, nbrs in self.edges.items():
            for j in nbrs:
                a, b = self.nodes[i], self.nodes[j]
                plt.plot([a[0], b[0]], [a[1], b[1]], "c--", alpha=0.3)

        # shortest path
        if self.path is not None:
            px = [p[0] for p in self.path]
            py = [p[1] for p in self.path]
            plt.plot(px, py, "r-", linewidth=3, label="Shortest Path")

        # start/goal
        plt.scatter(self.start[0], self.start[1], c="green", s=80, label="Start")
        plt.scatter(self.goal[0], self.goal[1], c="red", s=80, label="Goal")

        plt.axis("equal")
        plt.grid(True, alpha=0.3)
        plt.legend()

        if title is not None:
            plt.title(title)
        else:
            plt.title("Part 2: Visibility Graph Shortest Path")

        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()


# ===========================
# Dataset loading
# ===========================

def load_cases_from_file(fname):
    """
    Dataset format:

    # comment
    start_x start_y
    goal_x goal_y
    NUM_POLYGONS
    K1            # number of vertices in polygon 1
    x y
    ...
    K2
    ...

    Blank lines and comments (# ...) are ignored.
    Returns: list of dicts with keys: 'start', 'goal', 'polygons'.
    """
    # Force UTF-8 to avoid Windows GBK decode problems
    with open(fname, "r", encoding="utf-8") as f:
        raw_lines = f.readlines()

    # strip comments and blanks
    lines = []
    for line in raw_lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        lines.append(line)

    cases = []
    i = 0
    n = len(lines)
    while i < n:
        # start
        sx, sy = map(float, lines[i].split())
        i += 1
        # goal
        gx, gy = map(float, lines[i].split())
        i += 1
        # number of polygons
        num_polys = int(lines[i])
        i += 1

        polys = []
        for _ in range(num_polys):
            k = int(lines[i])  # number of vertices
            i += 1
            verts = []
            for _ in range(k):
                x, y = map(float, lines[i].split())
                i += 1
                verts.append([x, y])
            # ensure convex + ordered polygon using q1
            hull = findConvexHull(np.array(verts, float))
            polys.append(hull)

        cases.append(
            {
                "start": np.array([sx, sy], float),
                "goal": np.array([gx, gy], float),
                "polygons": polys,
            }
        )

    return cases


# ===========================
# Main
# ===========================

if __name__ == "__main__":
    if len(sys.argv) == 2:
        dataset_file = sys.argv[1]
        cases = load_cases_from_file(dataset_file)
        print(f"Loaded {len(cases)} cases from {dataset_file}")

        # create output folder for images
        out_dir = "results"
        os.makedirs(out_dir, exist_ok=True)

        for idx, case in enumerate(cases, start=1):
            print(f"\n=== Case {idx} ===")
            print("Start:", case["start"])
            print("Goal :", case["goal"])
            print("#Polygons:", len(case["polygons"]))

            vg = VisibilityGraph(case["polygons"], case["start"], case["goal"])

            img_path = os.path.join(out_dir, f"case_{idx}.png")

            if vg.path is None:
                print("No path found.")
            else:
                print("Path has", len(vg.path), "waypoints.")
                print("Saving image to:", img_path)

            vg.plot(
                title=f"Part 2 – Case {idx}",
                save_path=img_path,
            )

        print("\nAll figures saved in folder:", out_dir)

    else:
        # Fallback simple example if no dataset file is given
        print("Usage: python q2.py q2_dataset.txt")
        print("Running a simple built-in test instead...\n")

        polys = [
            findConvexHull([[1, 2], [4, 3], [4, 2]]),
            findConvexHull([[4, 8], [6, 7], [4, 4], [7, 6]]),
        ]
        start = [0, 0]
        goal = [10, 10]

        vg = VisibilityGraph(polys, start, goal)
        print("Shortest Path:", vg.path)
        vg.plot(title="Part 2 – Built-in Test")
