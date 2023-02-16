import numpy


def compute_RMSD(fixed, moving):
    if fixed.size == moving.size:
        rmsd = numpy.sqrt((((fixed - moving) ** 2) * 3).mean())

        return rmsd
