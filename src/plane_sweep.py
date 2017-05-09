import heapq
import logging
import math
import os

from bintrees import FastAVLTree

from src import gasv_input

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s|%(levelname)s: %(message)s')
logger.setLevel(logging.INFO)

INTERSECTION = 0
BEND = 1
END = 2
START = 3


def norm2(a):
    return (a[1], -a[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def left_turn(a, b, c):
    abx = b[0] - a[0]
    acx = c[0] - a[0]
    aby = b[1] - a[1]
    acy = c[1] - a[1]
    return abx * acy - acx * aby > 0


def intersect(s1, s2):
    return left_turn(s1[0], s1[1], s2[0]) != left_turn(s1[0], s1[1], s2[1]) and \
           left_turn(s2[0], s2[1], s1[0]) != left_turn(s2[0], s2[1], s1[1])


def intersection_point(s1, s2):
    dx = (s1[0][0] - s1[1][0], s2[0][0] - s2[1][0])
    dy = (s1[0][1] - s1[1][1], s2[0][1] - s2[1][1])

    div = det(dx, dy)
    assert (div != 0)

    d = (det(*s1), det(*s2))
    x = det(d, dx) / div
    y = det(d, dy) / div
    return x, y


class EventPointSchedule(object):
    def __init__(self, sweep_dir):
        self._sweep_dir = sweep_dir
        self.E = []

    def add_point(self, p):
        dst = dot(p, self._sweep_dir)
        heapq.heappush(self.E, (dst, p.type, p))

    def add_intersection(self, p):
        self.add_point(p)

    def has_next(self):
        return len(self.E) > 0

    def next(self):
        return heapq.heappop(self.E)

    def peek(self):
        return self.E[0]


def draw_polygons(name, reds, greens=[]):
    from PIL import Image, ImageDraw
    maxX = max(p[0] for p in sum(reds, []))
    minX = min(p[0] for p in sum(reds, []))
    maxY = max(p[1] for p in sum(reds, []))
    minY = min(p[1] for p in sum(reds, []))

    scale = 400.0 / (maxX - minX)
    img = Image.new("RGB", (400, math.ceil(scale * (maxY - minY))))
    drw = ImageDraw.Draw(img, 'RGBA')

    def _rescale(pl):
        return [((p[0] - minX) * scale, (maxY - p[1]) * scale) for p in pl]

    for i, p in enumerate(reds):
        drw.polygon(_rescale(p), (255, 0, i * 20, int(55 + 200 / len(reds))),
                    outline=(0, 0, 255, 125))
    for i, p in enumerate(greens):
        drw.polygon(_rescale(p), (0, 255, 0, int(55 + 200 / len(greens))))
    del drw

    img.save('%s.png' % (name), 'PNG')


class SweepLineStatus(object):
    def __init__(self, sweep_dir):
        self._tree = FastAVLTree()
        self._sweep_dir = sweep_dir

    def set_position(self, pos):
        self._pos = pos
        L0 = pos * self._sweep_dir
        L1 = norm2(self._sweep_dir)
        self._line = (L0, L1)

    def add_edge(self, s, d):
        self._tree.insert(dot(self._line, s), d)

    def find(self, p):
        return (self._tree.succ_item(dot(self._line, p)),
                self._tree.prev_item(dot(self._line, p)))

    def remove(self, p):
        self._tree.remove(p)


class PlaneSweep(object):
    def __init__(self, sweep_dir, min_intersections=4, output_directory=None):
        self._sweep_dir = sweep_dir
        self._schedule = EventPointSchedule(sweep_dir)
        self._intersections = []
        self._min_intersections = min_intersections
        self._output_directory = output_directory

    def add_polygon(self, polygon: gasv_input.Polygon):
        for i, e in enumerate(polygon.edges):
            pA = polygon.edges[i - 1].p1
            p = e.p1
            pB = e.p2

            p_d = dot(self._sweep_dir, p)
            pA_d = dot(self._sweep_dir, pA)
            pB_d = dot(self._sweep_dir, pB)

            p.outgoing = []
            p.incoming = []

            if pA_d > p_d:
                p.outgoing.append(pA)
            else:
                p.incoming.append(pA)
            if pB_d >= p_d:
                p.outgoing.append(pB)
            else:
                p.incoming.append(pB)

            if len(p.outgoing) > 0 and len(p.incoming) > 0:
                p.type = BEND
            elif len(p.outgoing) > 0:
                p.type = START
            else:
                p.type = END
            p.poly = polygon
            p.polyset = {polygon}

            self._schedule.add_point(p)

    def _handle_intersection(self, P, L, A, B, As, Bs, d):
        for s in P.outgoing:
            t = (P, s)
            for idx, l in enumerate(L):
                if intersect(t, l):
                    P_x = intersection_point(t, l)
                    P_i = gasv_input.Point((P_x[0], P_x[1]))
                    P_i.type = INTERSECTION
                    P_i.poly = None
                    P_i.polyset = P.polyset | l[0].polyset
                    P_i.incoming = [P, l[0]]
                    P_i.outgoing = [s, l[1]]

                    assert d <= dot(self._sweep_dir, P_i)
                    self._schedule.add_intersection(P_i)
                    assert (P.poly != l[0].poly)

    def find_intersections(self):
        # L = SweepLineStatus(self._sweep_dir)
        L = []
        A = []
        B = []
        Ra = []
        Rb = []
        while self._schedule.has_next():
            d, _, P = self._schedule.next()

            if P.type == START:
                self._handle_intersection(P, L, A, B, [], [], d)

                t = left_turn(P, P.outgoing[0], P.outgoing[1])

                for i in P.outgoing:
                    assert d <= dot(self._sweep_dir, i)
                    L.append((P, i))
                    A.append([P])
                    B.append([P])
                    if t:
                        Ra.append(P.polyset)
                        Rb.append(set())
                    else:
                        Ra.append(set())
                        Rb.append(P.polyset)
                    t = not t
            elif P.type == BEND:
                i = P.incoming[0]
                idx = L.index((i, P))
                L.pop(idx)
                prev_a = A.pop(idx)
                prev_b = B.pop(idx)
                r_a = Ra.pop(idx)
                r_b = Rb.pop(idx)
                self._handle_intersection(P, L, A, B, prev_a, prev_b, d)
                for i in P.outgoing:
                    assert d <= dot(self._sweep_dir, i)
                    L.append((P, i))
                    A.append(prev_a + [P])
                    B.append([P] + prev_b)
                    Ra.append(r_a)
                    Rb.append(r_b)
            elif P.type == END:
                s = (P.incoming[0], P)
                t = (P.incoming[1], P)

                if not left_turn(P, s[0], t[0]):
                    s, t = t, s
                idxS = L.index(s)
                idxT = L.index(t)

                concatx = A[idxT] + [P] + B[idxS]
                # concaty = A[idxS] + [P] + B[idxT]

                for x in (concatx,):
                    parents = set.union(Rb[idxS], Ra[idxT])
                    if len(parents) >= self._min_intersections:
                        self.report_intersection(parents, x)
                if idxS < idxT:
                    idxS, idxT = idxT, idxS
                L.pop(idxS)
                A.pop(idxS)
                B.pop(idxS)
                Ra.pop(idxS)
                Rb.pop(idxS)
                L.pop(idxT)
                A.pop(idxT)
                B.pop(idxT)
                Ra.pop(idxT)
                Rb.pop(idxT)
            elif P.type == INTERSECTION:
                s = (P.incoming[0], P.outgoing[0])
                t = (P.incoming[1], P.outgoing[1])

                if not left_turn(P, s[0], t[0]):
                    s, t = t, s
                idxS = L.index(s)
                idxT = L.index(t)

                x = A[idxT] + [P] + B[idxS]
                parents = set.union(Rb[idxS], Ra[idxT])
                if len(parents) >= self._min_intersections:
                    self.report_intersection(parents, x)
                A_s = A[idxS]
                A[idxS] = [P]
                B[idxS] = [P] + B[idxT]
                A[idxT] = A_s + [P]
                B[idxT] = [P]

                # Update region tracking status
                r_a = Ra[idxS]
                r_b = Rb[idxS]
                Ra[idxS] = (Ra[idxS] | (Rb[idxT] - Ra[idxT])) - (Ra[idxT] - Rb[idxT])
                Rb[idxS] = (Rb[idxS] | (Rb[idxT] - Ra[idxT])) - (Ra[idxT] - Rb[idxT])
                Ra[idxT] = (Ra[idxT] | (r_a - r_b)) - (r_b - r_a)
                Rb[idxT] = (Rb[idxT] | (r_a - r_b)) - (r_b - r_a)

    def report_intersection(self, parents, points):
        parent_names = sorted([p.breakpoint.name for p in parents])
        logger.info("Intersection: region=%s, reads=%s", points, " ".join(parent_names))
        self._intersections.append((parents, points))
        if self._output_directory:
            os.makedirs(self._output_directory, exist_ok=True)
            draw_polygons(
                os.path.join(self._output_directory, "#".join(parent_names)),
                [p.points for p in parents], [points])

    def intersections(self):
        return self._intersections


def find_maximal_intersections(intersections):
    # for R, p in intersections:
    return intersections


def find_intersections(breakpoints, sweep_dir, min_intersections, output_directory):
    ps = PlaneSweep(sweep_dir, min_intersections, output_directory)
    for bp in breakpoints:
        ps.add_polygon(bp.polygon)

    ps.find_intersections()
    return find_maximal_intersections(ps.intersections())
