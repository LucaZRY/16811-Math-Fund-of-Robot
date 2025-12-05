
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

from q1 import findConvexHull
from q2 import VisibilityGraph



def minkowski_sum_convex(obstacle, robot):
    """
    Compute convex Minkowski sum: obstacle ⊕ (-robot)
    i.e. all points o - r, o in obstacle, r in robot.
    Brute force all pairwise sums, then convex hull.
    obstacle, robot: (N x 2) numpy arrays (convex).
    Returns: (M x 2) convex hull vertices (C-obstacle).
    """
    pts = []
    for o in obstacle:
        for r in robot:
            pts.append(o - r)
    pts = np.array(pts, dtype=float)
    hull = findConvexHull(pts)
    return hull


def build_cspace_obstacles(workspace_obstacles, robot_poly):
    """
    workspace_obstacles: list of convex polygons (workspace).
    robot_poly: convex polygon in its own local frame (reference at (0,0)).
    Returns list of convex polygons in configuration space.
    """
    cspace_polys = []
    for obs in workspace_obstacles:
        c_obs = minkowski_sum_convex(obs, robot_poly)
        cspace_polys.append(c_obs)
    return cspace_polys


def translate_polygon(poly, t):
    """
    Translate polygon poly (N x 2) by vector t (2,).
    """
    poly = np.asarray(poly, float)
    t = np.asarray(t, float)
    return poly + t


def plot_workspace(robot_poly, workspace_obstacles, path, save_path=None, title=None):
    """
    Plot original workspace:
    - workspace obstacles
    - robot path as line (reference point)
    - robot polygon at start and goal configurations.
    """
    plt.figure()

    for P in workspace_obstacles:
        P = np.asarray(P)
        xs = np.append(P[:, 0], P[0, 0])
        ys = np.append(P[:, 1], P[0, 1])
        plt.plot(xs, ys, "k-")

    if path is not None and len(path) >= 2:
        px = [p[0] for p in path]
        py = [p[1] for p in path]
        plt.plot(px, py, "r-", linewidth=2, label="Robot reference path")

        plt.scatter(px[0], py[0], c="green", s=80, label="Start config")
        plt.scatter(px[-1], py[-1], c="blue", s=80, label="Goal config")

        R_start = translate_polygon(robot_poly, path[0])
        R_goal = translate_polygon(robot_poly, path[-1])

        for poly, color, lab in [(R_start, "g", "Robot at start"),
                                 (R_goal, "b", "Robot at goal")]:
            xs = np.append(poly[:, 0], poly[0, 0])
            ys = np.append(poly[:, 1], poly[0, 1])
            plt.plot(xs, ys, color + "-", linewidth=2, label=lab)

    plt.axis("equal")
    plt.grid(True, alpha=0.3)
    plt.legend()

    if title is not None:
        plt.title(title)
    else:
        plt.title("Part 3: Workspace (Robot + Obstacles)")

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def load_q3_cases(fname):
    """
    q3 dataset format (q3_dataset.txt):

    ROBOT
    K_robot
    rx1 ry1
    ...
    rxK ryK

    # Case 1
    start_x start_y
    goal_x goal_y
    NUM_OBSTACLES
    K1
    ox11 oy11
    ...
    K2
    ...

    # Case 2
    start_x start_y
    ...

    Lines beginning with '#' and blank lines are ignored.

    Returns:
      robot_poly: (R x 2) numpy array (convex hull of robot vertices)
      cases: list of dicts
        {
          "start": np.array([sx, sy]),
          "goal" : np.array([gx, gy]),
          "obstacles": [poly1, poly2, ...]  # workspace polys (convex)
        }
    """
    with open(fname, "r", encoding="utf-8") as f:
        raw_lines = f.readlines()

    lines = []
    for line in raw_lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        lines.append(line)

    if not lines or lines[0].upper() != "ROBOT":
        raise ValueError("Dataset must start with line 'ROBOT'")

    i = 1
    n = len(lines)

    if i >= n:
        raise ValueError("Missing robot vertex count after 'ROBOT'")
    k_robot = int(lines[i])
    i += 1
    if i + k_robot > n:
        raise ValueError("Not enough lines for robot vertices")

    robot_verts = []
    for _ in range(k_robot):
        x, y = map(float, lines[i].split())
        i += 1
        robot_verts.append([x, y])
    robot_verts = np.array(robot_verts, float)
    robot_poly = findConvexHull(robot_verts)  

    cases = []
    while i < n:
        sx, sy = map(float, lines[i].split())
        i += 1
        gx, gy = map(float, lines[i].split())
        i += 1
        num_obs = int(lines[i])
        i += 1

        obstacles = []
        for _ in range(num_obs):
            k = int(lines[i]) 
            i += 1
            verts = []
            for _ in range(k):
                x, y = map(float, lines[i].split())
                i += 1
                verts.append([x, y])
            verts = np.array(verts, float)
            obs_poly = findConvexHull(verts)  
            obstacles.append(obs_poly)

        cases.append(
            {
                "start": np.array([sx, sy], float),
                "goal": np.array([gx, gy], float),
                "obstacles": obstacles,
            }
        )

    return robot_poly, cases

if __name__ == "__main__":
    if len(sys.argv) == 2:
        dataset_file = sys.argv[1]
        robot_poly, cases = load_q3_cases(dataset_file)
        print(f"Loaded robot + {len(cases)} cases from {dataset_file}")

        out_dir_cspace = "results_q3_cspace"
        out_dir_ws = "results_q3_workspace"
        os.makedirs(out_dir_cspace, exist_ok=True)
        os.makedirs(out_dir_ws, exist_ok=True)

        for idx, case in enumerate(cases, start=1):
            print(f"\n=== Case {idx} ===")
            print("Start config:", case["start"])
            print("Goal  config:", case["goal"])
            print("#Workspace obstacles:", len(case["obstacles"]))

            cspace_polys = build_cspace_obstacles(case["obstacles"], robot_poly)

            vg = VisibilityGraph(cspace_polys, case["start"], case["goal"])
            path = vg.path  

            if path is None:
                print("No collision-free path found in configuration space.")
                c_img = os.path.join(out_dir_cspace, f"case_{idx}_cspace.png")
                vg.plot(title=f"Part 3 – Case {idx} (C-space, no path)", save_path=c_img)
                continue

            print("Found path with", len(path), "waypoints.")

            c_img = os.path.join(out_dir_cspace, f"case_{idx}_cspace.png")
            vg.plot(title=f"Part 3 – Case {idx} (Configuration Space)", save_path=c_img)

            ws_img = os.path.join(out_dir_ws, f"case_{idx}_workspace.png")
            plot_workspace(
                robot_poly,
                case["obstacles"],
                path,
                save_path=ws_img,
                title=f"Part 3 – Case {idx} (Workspace)",
            )

        print("\nC-space figures saved in:", out_dir_cspace)
        print("Workspace figures saved in:", out_dir_ws)

    else:
        print("Usage: python q3.py q3_dataset.txt")
        print("Running a simple built-in demo instead...\n")

        robot_poly = findConvexHull(
            np.array(
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [1.0, 0.5],
                    [0.0, 0.5],
                ],
                float,
            )
        )

        obstacles = [
            findConvexHull(
                np.array(
                    [
                        [3.0, 2.0],
                        [5.0, 2.0],
                        [5.0, 4.0],
                        [3.0, 4.0],
                    ],
                    float,
                )
            )
        ]

        start = np.array([0.0, 0.5])
        goal = np.array([8.0, 0.5])

        cspace_polys = build_cspace_obstacles(obstacles, robot_poly)
        vg = VisibilityGraph(cspace_polys, start, goal)
        path = vg.path

        if path is None:
            print("No path found in demo.")
        else:
            print("Demo path waypoints:")
            for p in path:
                print(p)

        vg.plot(title="Part 3 – Demo (Configuration Space)")
        plot_workspace(robot_poly, obstacles, path, title="Part 3 – Demo (Workspace)")
