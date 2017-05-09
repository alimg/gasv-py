import argparse
import logging
import os

from src import gasv_input
from src import plane_sweep

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s|%(levelname)s: %(message)s')
logger.setLevel(logging.INFO)


def read_gasv_files(breakpoints):
    pass


def run():
    parser = argparse.ArgumentParser(description='GASV-PY')
    parser.add_argument('input',
                    help='load GASV breakpoints file (e.g. Example.bam.gasv.in)')
    parser.add_argument('--min-cluster-size', type=int, default=4,
                        help='set minimum number of intersections per cluster')
    args = parser.parse_args()

    with open(args.input) as input_file:
        for line in input_file:
            file_name, read_type, Lmin, Lmax = line.rstrip().split("\t")
            if read_type != "PR":
                logger.error("unsupported file type '%s'. skipping" % read_type)
                continue
            logger.info("processing '%s'" % file_name)
            breakpoints = list(gasv_input.read_gasv_file(file_name, int(Lmin), int(Lmax)))
            results = plane_sweep.find_intersections(
                breakpoints, (1, 1.0001), args.min_cluster_size,
                output_directory=os.path.join("out", os.path.basename(file_name)))
            logger.info("found %d regions", len(results))


if __name__ == "__main__":
    run()
