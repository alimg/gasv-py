
class Point(tuple):
    pass


class Edge(object):

    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2


class Polygon(object):
    def __init__(self, breakpoint, points):
        self.breakpoint = breakpoint
        self.points = [Point(p) for p in points]
        self.edges = self._edges()

    def _edges(self):
        edges = []
        for i, p in enumerate(self.points):
            edges.append(Edge(self.points[i-1], p))
        return edges


class BreakPoint(object):
    def __init__(self, name, xA, xB, xS, yA, yB, yS, lmin, lmax):
        self.name = name
        self.x_start = int(xS + xA)
        self.x_end = int(xS + xB)
        self.xS = xS
        self.y_start = int(yS + yA)
        self.y_end = int(yS + yB)
        self.yS = yS
        self.lmin = lmin
        self.lmax = lmax

        self.x = min(self.x_start, self.x_end)
        self.y = min(self.y_start, self.y_end)

        self.polygon = self._polygon()

    def _polygon(self):
        dirX = 1 if self.xS == '+' else -1
        dirY = 1 if self.yS == '+' else -1
        x = abs(self.x)
        y = abs(self.y)

        return Polygon(self, [
            (x + dirX * self.lmin, y),
            (x + dirX * self.lmax, y),
            (x, y + dirY * self.lmax),
            (x, y + dirY * self.lmin)])


def read_gasv_file(filename, lmin, lmax):
    with open(filename, "r") as rf:
        for line in rf:
            row = line.strip().split("\t")
            yield BreakPoint(row[0], row[2], row[3], row[4],
                             row[6], row[7], row[8], lmin, lmax)
