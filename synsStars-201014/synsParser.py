
import argparse
import clDatabase

def getParser():
# {{{
    """Returns parser for the synStars code."""
    
    # setup parser
    parser = argparse.ArgumentParser(
    parents = [clDatabase.getParser()], 
    add_help=False)
    parser.add_argument("--mtot", type=float, help="Total mass of the cluster in MSun.")
    parser.add_argument("--tmax", default=40., type=float, help="Maximum time in Myr. Default is 40 Myr.")
    parser.add_argument("--Nmbin", default=5000, type=int, help="Number of mass bins over ranges of imfbin. Mass of the cluster, mtot, Default is 5000.")
    parser.add_argument("--Nt", default=500, type=int, help="Number points of time evolution (0 to tmax). Default is 100.")
    parser.add_argument("--imfbin", default='0.01,0.08,0.5,120', type=str, help="IMF mass bins in [MSun], coma separated list. \
    Mass of the cluster, mtot, is distributed over these all bins. However, track-grid is created only between provided tracks. \
    For mass ranges where no tracks are available, mass and energy deposition rates are assumed to be zero. Default is 0.01,0.08,0.5,120.")
    parser.add_argument("--imftype", type=str, default="ppl", help="Type of the IMF \
    salpeter is a single power-law with slope -2.35; minimum and maximum mass can be controlled by -imfbin. \
    kroupa is piecewise power-law with slope -0.3 below 0.08 MSun, -1.3 between 0.08 and 0.5 MSun and -2.3 above 0.5 MSun; \
    minimum and maximum mass can be controlled by -imfbin. \
    ppl is a general piecewise power-law controlled with arguments -imfbin and -alpha. \
    maschberger is a continuous heavy-tailed log-normal distribution with the slopes of the tails close to -0.48 \
    below 0.043 MSun and -2.3 above 0.93 MSun. The function can be controlled by arguments \
    -alpha and -imfbin. \
    default value is ppl with the alpha values equal to Kroupa.")
    parser.add_argument("--alpha", type=str, default="-0.3,-1.3,-2.3", help="Indeces of piecewise power-law IMF, coma separated list. \
    Shortcuts kroupa and salpeter can also be used. \
    Default value is kroupa.")
    parser.add_argument("--ret1", default=1.0, type=float, help="Retention fraction of the 1st stellar generation. Default is 1.")
    parser.add_argument("--imffile", default="", type=str, help="File template for writing imf evolution. Default is '', i.e. not written.")
    parser.add_argument("--trdir", default="./tracks/", type=str, help="Directory with track files (*.track.dat). Default is './tracks/'.")
    parser.add_argument("--trfile", default=".*\.track\.dat$", type=str, help="Track file name template regexp. Default is '.*\.track\.dat'.")
    parser.add_argument("--gridfile", default="", type=str, help="File with the interpolated fine grid of tracks. Default is '', i.e. file is not written.")
    parser.add_argument("--intType_tM0", default="log", type=str, help="Type of interpolation (log or lin) of the time as a function of ZAMS mass M0. \
    By default, this interpolation is done in the logarithmic scale, i.e. log(t) is interpolated as a function of log(M0)")
    parser.add_argument("--iord_tM0", default=2, type=int, help="Order of the interpolation (spline) for the time(M_ZAMS) function. Default is 2.")
    parser.add_argument("--iord_QM0", default=1, type=int, help="Order of the interpolation (spline) for all (other than time) quantities as functions of M_ZAMS. Default is 1.")

    return parser
# }}}

