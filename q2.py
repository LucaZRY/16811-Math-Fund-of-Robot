import sys
import itertools
import numpy as np
from q1 import findConvexHull
from matplotlib import pyplot as plt


class VisibiltyGraph:
    """
    Visibility graph + Dijkstra shortest path for:
      - Part 2: point robot (use polygons as given)
      - Part 3: polygonal robot (call getMinkowskiSum(robot) first to build C-space)

    Usage pattern:
        vg = VisibiltyGraph(points, start, end)
        # for part 2 (point robot) just:
        vg.findShortestPath()
        vg.plotPolygonsAndPaths()

        # for part 3 (polygon robot):
        vg.getMinkowskiSum(robot_vertices)
        vg.findShortestPath()
        vg.plotPolygonsAndPaths(robot_vertices, isRobot=True)
    """

    def __init__(self, points, start, end):
        self.polygons = []
        self.original_polygons = []

        # Assign the start and end for the robot to move
        self.start = np.array(start, dtype=float)
        self.end = np.array(end, dtype=float)
        self.shortestPath = None

        # Form the initial set of polygons
        self.formPolygonsFromPoints(points)
        # Save for plotting
        self.original_polygons = list(self.polygons)
        # Compute visibility graph (with overlap handling)
        self.checkAndCompute()

    # ------------------------------------------------------------------
    # Polygon and overlap handling
    # ------------------------------------------------------------------

    def formPolygonsFromPoints(self, points):
        # form convex polygons and assign them indices
        for i, set_points in enumerate(points):
            set_points = np.array(set_points, np.float32)
            chull = findConvexHull(set_points.copy())
            poly_elem = {"convex_hull": chull, "points": set_points, 'index': i}
            self.polygons.append(poly_elem)

    def checkAndCompute(self):
        # ensure valid start and end
        self.checkForBaseConditions()
        # check for overlapping polygons
        isOverlap = self.checkOverlappingPolygons()
        if isOverlap:
            self.recreateOverlappingPolygons()
        # compute visibility graph
        self.formVisibiltyGraph()

    def findSetOfTwoPoints(self, points):
        # return 2 points in order starting from index 0
        indices = [i for i in range(points.shape[0])] + [0]
        for i in range(len(indices) - 1):
            val = indices[i:i + 2]
            if len(val) == 2:
                yield val

    def isPointInGivenPolygon(self, point, polygon):
        # ray-casting point-in-polygon test, for convex hull
        x, y = point
        inside = False
        hull_points = polygon['convex_hull']

        n = len(hull_points)
        p1x, p1y = hull_points[0][0], hull_points[0][1]
        for i in range(0, n + 1):
            p2x, p2y = hull_points[i % n][0], hull_points[i % n][1]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y

        return inside

    def arePolygonsIntersecting(self, poly1, poly2):
        # check if given two polys intersect by checking hull vertices
        for point in poly1['convex_hull']:
            if self.isPointInGivenPolygon(point, poly2):
                return True
        return False

    def checkOverlappingPolygons(self):
        # check if any two polys are overlapping
        for (poly1, poly2) in itertools.combinations(self.polygons, 2):
            if self.arePolygonsIntersecting(poly1, poly2):
                return True
        return False

    def checkForBaseConditions(self):
        # start and end must not lie inside any polygon
        for polygon in self.polygons:
            if self.isPointInGivenPolygon(self.start, polygon):
                print('No Valid Path exists! (start inside obstacle)')
                sys.exit(1)
            if self.isPointInGivenPolygon(self.end, polygon):
                print('No Valid Path exists! (end inside obstacle)')
                sys.exit(1)

    def recreateOverlappingPolygons(self):
        # Merge overlapping polygons into single convex obstacles
        new_polygons = []

        intersecting_polys = {self.polygons[i]['index']: set()
                              for i in range(len(self.polygons))}
        for (poly1, poly2) in itertools.combinations(self.polygons, 2):
            if self.arePolygonsIntersecting(poly1, poly2):
                intersecting_polys[poly1['index']].add(poly2['index'])
                intersecting_polys[poly2['index']].add(poly1['index'])

        for k, v in intersecting_polys.items():
            if len(v) == 0:
                new_polygons.append(self.polygons[k])
            else:
                _point = self.fetchHullPoints(k)
                new_index = [k]
                for _v in list(v):
                    _point = np.vstack((_point, self.fetchHullPoints(_v)))
                    new_index.append(_v)

                new_index.sort(reverse=True)
                new_index = [str(a) for a in new_index]
                new_index = int(''.join(new_index))
                hull = findConvexHull(_point)
                new_polygons.append(
                    {"convex_hull": hull, "points": _point, 'index': new_index}
                )

        self.polygons = []
        q_dict = {}
        for pol in new_polygons:
            if pol['index'] in q_dict:
                continue
            q_dict[pol['index']] = 1
            self.polygons.append(pol)

        overlap = self.checkOverlappingPolygons()
        if overlap:
            self.formVisibiltyGraph()
            self.recreateOverlappingPolygons()

    def fetchHullPoints(self, index):
        for i in self.polygons:
            if i['index'] == index:
                return i['convex_hull']

    # ------------------------------------------------------------------
    # Visibility graph construction
    # ------------------------------------------------------------------

    def areConnectedNeighbors(self, v1, v2, e1, e2):
        # check if the vertices are connected to neighboring edges
        if (v2 == e2).all() and not (v1 == e1).all():
            return True
        if (v2 == e1).all() and not (v1 == e2).all():
            return True
        if (v1 == e2).all() and not (v2 == e2).all():
            return True
        if (v1 == e1).all() and not (v2 == e2).all():
            return True
        return False

    def getPowerofPoint(self, line, point):
        # side test of a point wrt a line (e1,e2)
        e1, e2 = line[0], line[1]
        m = e1[1] - e2[1]
        if (e1[0] - e2[0]) == 0:
            return point[1] - e1[0]

        m /= (e1[0] - e2[0])
        c = e1[1] - m * e1[0]
        side = m * point[0] - point[1] + c
        side = side / np.sqrt(1 + m ** 2)
        if np.abs(side) < 1e-4:
            side = 0
        sign = np.sign(side)
        return sign

    def findPolygonalConnectionsinVisibiltyGraph(self):
        # add polygon adjacency edges (neighbors on each hull) into visibility graph
        for k, _ in self.visibility_graph.items():
            temp_k = tuple(k)
            polygon_id = int(self.vertex_map[temp_k])
            try:
                poly_pts = self.polygons[polygon_id]['convex_hull']
            except Exception:
                print('No valid path exist as no visibility vertex is found!')
                sys.exit(1)
            lngth = poly_pts.shape[0]

            for i in range(poly_pts.shape[0]):
                if (poly_pts[i] == k).all():
                    self.visibility_graph[temp_k].add(tuple(poly_pts[(i - 1) % lngth]))
                    self.visibility_graph[temp_k].add(tuple(poly_pts[(i + 1) % lngth]))
                    break

    def formVisibiltyGraph(self):
        if len(self.polygons) == 0:
            print('No obstacles. Shortest path is a straight line from {0} to {1}'.format(
                self.start, self.end))
            sys.exit(1)

        # load the points per poly (with a polygon index)
        temp_graph = []
        for polynum, poly in enumerate(self.polygons):
            for point in poly['convex_hull']:
                index = np.array([polynum])
                a = np.hstack((point, index))
                temp_graph.append(a)

        # add the source and end vertex as they dont exist in any poly
        temp_graph.append(np.array([self.start[0], self.start[1], -1]))
        temp_graph.append(np.array([self.end[0], self.end[1], -2]))
        temp_graph = np.array(temp_graph)

        # create a map of vertices -> polygon index (or -1/-2)
        self.vertex_map = {}
        for _ in self.polygons:
            for vertex in temp_graph:
                self.vertex_map[tuple(vertex[:-1])] = vertex[-1]

        self.inverse_vertex_map = {v: k for k, v in self.vertex_map.items()}

        visibility_graph = {}
        for g in temp_graph:
            # for each vertex, add a possible visible graph
            visibility_graph[tuple(g[:-1])] = set()

        # compute polygon edges
        edges = []
        for poly in self.polygons:
            for _, (i1, i2) in enumerate(self.findSetOfTwoPoints(poly['convex_hull'])):
                p1, p2 = poly['convex_hull'][i1], poly['convex_hull'][i2]
                if (p1 == p2).all():
                    continue
                edges.append((p1, p2))

        edges = np.array(edges)

        # for each pair of vertices, check visibility against all edges
        for _, vertex1 in enumerate(temp_graph):
            id1 = vertex1[-1]
            vertex1 = vertex1[:-1]
            for _, vertex2 in enumerate(temp_graph):
                id2 = vertex2[-1]
                vertex2 = vertex2[:-1]
                if id1 == id2:
                    break
                if (vertex1 == vertex2).all():
                    continue

                blocked = False
                for _, (edge1, edge2) in enumerate(edges):
                    # skip degenerate and duplicate edges
                    if (edge1 == edge2).all():
                        continue
                    if ((edge1 == vertex1).all() and (edge2 == vertex2).all()) or \
                       ((edge1 == vertex2).all() and (edge2 == vertex1).all()):
                        continue
                    if self.areConnectedNeighbors(vertex1, vertex2, edge1, edge2):
                        continue

                    x1, y1 = vertex1
                    x2, y2 = vertex2
                    x3, y3 = edge1
                    x4, y4 = edge2
                    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
                    if denom != 0:
                        px = ((x1 * y2 - y1 * x2) * (x3 - x4) -
                              (x1 - x2) * (x3 * y4 - y3 * x4)) / denom
                        py = ((x1 * y2 - y1 * x2) * (y3 - y4) -
                              (y1 - y2) * (x3 * y4 - y3 * x4)) / denom

                        list_vertex = [vertex1, vertex2]
                        if self.getPowerofPoint(list_vertex, edge1) == self.getPowerofPoint(list_vertex, edge2):
                            # line is on same side w.r.t edge; no crossing
                            continue

                        if (px >= min(vertex1[0], vertex2[0]) and
                            px <= max(vertex1[0], vertex2[0]) and
                            py >= min(vertex1[1], vertex2[1]) and
                            py <= max(vertex1[1], vertex2[1])):
                            # segment intersects this edge in its interior → blocked
                            blocked = True
                            break
                if not blocked:
                    visibility_graph[tuple(vertex1)].add(tuple(vertex2))
                    visibility_graph[tuple(vertex2)].add(tuple(vertex1))

        self.visibility_graph = visibility_graph
        self.findPolygonalConnectionsinVisibiltyGraph()

    # ------------------------------------------------------------------
    # Dijkstra / shortest path
    # ------------------------------------------------------------------

    def createAdjacencyMatrix(self):
        keys = list(self.visibility_graph.keys())

        adjMatrix = np.full(shape=(len(keys), len(keys)), fill_value=np.inf)
        for k, v in self.visibility_graph.items():
            key_index = keys.index(k)
            for vertices in v:
                v_index = keys.index(vertices)
                weight = np.linalg.norm(np.array(k) - np.array(vertices))
                adjMatrix[key_index, v_index] = weight
                adjMatrix[v_index, key_index] = weight

        self.keys = keys

        source_index = keys.index((self.start[0], self.start[1]))
        target_index = keys.index((self.end[0], self.end[1]))

        # explicitly keep them as unknown in the matrix; path is determined by Dijkstra
        adjMatrix[source_index, target_index] = np.inf
        adjMatrix[target_index, source_index] = np.inf

        return adjMatrix, source_index, target_index

    def findShortestPath(self):
        final_matrix, source, target = self.createAdjacencyMatrix()

        num_vertices = len(self.keys)
        dist = [np.inf for _ in range(num_vertices)]
        prev = [None for _ in range(num_vertices)]

        dist[source] = 0
        Q = [i for i in range(num_vertices)]

        def findMininQueue(Q, dist):
            min_dist = np.inf
            min_index = None
            for q in Q:
                d = dist[q]
                if d < min_dist:
                    min_dist = d
                    min_index = q
            return min_index

        visited = set()
        while len(Q):
            u = findMininQueue(Q, dist)
            visited.add(u)
            Q.remove(u)
            neighbors = np.where(final_matrix[u] != np.inf)[0]
            for v in neighbors:
                if v not in visited:
                    alt = dist[u] + final_matrix[u, v]
                    if alt < dist[v]:
                        dist[v] = alt
                        prev[v] = u

        # reconstruct path from target to source and reverse it
        u = target
        path = [u]
        while u != source:
            u = prev[u]
            path.append(u)
        path.reverse()

        actual_path = [self.keys[i] for i in path]
        self.shortestPath = actual_path

    # ------------------------------------------------------------------
    # Minkowski sum for Part 3
    # ------------------------------------------------------------------

    def getMinkowskiSum(self, robot):
        """
        Expand obstacles by Minkowski sum with a convex robot:
            O_i ⊕ (−R)

        This turns the polygonal robot problem into a point-robot
        problem in configuration space.
        """
        crobot = findConvexHull(robot)
        robot = -1 * crobot

        # configuration-space obstacles
        curr_polys = []
        for polygon in self.polygons:
            curr_poly = []
            for poly_pt in polygon['convex_hull']:
                for robot_pt in robot:
                    curr_poly.append(poly_pt + robot_pt)
            curr_polys.append(curr_poly)

        final_polys = []
        for curr_index, poly_points in enumerate(curr_polys):
            chull = findConvexHull(poly_points)
            poly_elem = {"convex_hull": chull,
                         "points": poly_points,
                         'index': curr_index}
            final_polys.append(poly_elem)

        self.polygons = final_polys
        # recompute visibility graph / checks in C-space
        self.checkAndCompute()

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def plotPolygonsAndPaths(self, robot=None, isRobot=False):
        # plot original workspace polygons (not expanded C-space ones)
        for poly in self.original_polygons:
            poly_hull = poly['convex_hull']
            hull_x = poly_hull[:, 0]
            hull_y = poly_hull[:, 1]

            hull_x = np.append(hull_x, hull_x[0])
            hull_y = np.append(hull_y, hull_y[0])
            plt.plot(hull_x, hull_y, "b-")

        if self.shortestPath:
            spath_x = [float(p[0]) for p in self.shortestPath]
            spath_y = [float(p[1]) for p in self.shortestPath]

            for i in range(len(self.shortestPath) - 1):
                plt.plot([spath_x[i], spath_x[i + 1]],
                         [spath_y[i], spath_y[i + 1]], 'k-')

        plt.plot(self.start[0], self.start[1], 'ro')
        plt.plot(self.end[0], self.end[1], 'go')

        if isRobot:
            if robot is None:
                print('Pass the robot for plotting')
                sys.exit(1)
            else:
                plot_robo = np.array(robot)
                hull_x = plot_robo[:, 0]
                hull_y = plot_robo[:, 1]
                hull_x = np.append(hull_x, hull_x[0])
                hull_y = np.append(hull_y, hull_y[0])
                plt.plot(hull_x, hull_y, "g-")

        plt.grid(True)
        plt.axis("equal")
        plt.show()


if __name__ == "__main__":
    # Example test (point robot, Part 2 style)
    points = [
        [[0.0, 1.0], [4.0, 2.0], [4.0, 0.0]],
        [[4.0, 8.0], [6.0, 7.0], [4.0, 5.0], [7.0, 6.0]],
        [[6.0, 8.0], [9.84, 8.87], [7.16, 11.76]],
        [[6.0, 10.0], [8.0, 10.0], [5.0, 11.74],
         [6.0, 7.46], [8.0, 13.46], [9.0, 11.73]],
    ]

    start = np.array([0, 0], dtype=float)
    end = np.array([8, 8], dtype=float)

    vg = VisibiltyGraph(points, start, end)
    vg.findShortestPath()
    print("Shortest path:", vg.shortestPath)
    vg.plotPolygonsAndPaths()
