from src import gasv_input
from src import plane_sweep


def run():
    Lmin = 153
    Lmax = 247
    breakpoints = list(gasv_input.read_gasv_file("Example.bam_all.deletion", Lmin, Lmax))
    # breakpoints = list(gasv_input.read_gasv_file("Example.bam_all.inversion", Lmin, Lmax))
    print(breakpoints[0].polygon.points)

    plane_sweep.find_intersections(breakpoints, (1, 1.0001))


if __name__ == "__main__":
    run()
