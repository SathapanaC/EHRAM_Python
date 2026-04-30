import sys
import argparse
import numpy as np
from datetime import datetime

class ErhamState:
    """
    Class to hold the shared state equivalent to the Fortran COMMON blocks.
    """
    def __init__(self):
        # Arrays based on COMMON blocks
        self.d = np.zeros((2, 243, 243), dtype=np.float64)
        self.d1 = np.zeros((243, 243), dtype=np.float64)
        self.d2 = np.zeros((243, 243), dtype=np.float64)
        self.phi1 = np.zeros(970, dtype=np.float64)
        self.phi2 = np.zeros(485, dtype=np.float64)
        self.c = np.zeros((243, 243), dtype=np.float64)
        self.u = np.zeros((2, 243, 243), dtype=np.float64)
        self.e = np.zeros(243, dtype=np.float64)
        self.ej = np.zeros(243, dtype=np.float64)
        self.h = np.zeros((2, 243, 243), dtype=np.float64)
        self.ew = np.zeros((243, 243), dtype=np.float64)
        self.b = np.zeros((2, 243, 243), dtype=np.float64)
        self.ev = np.zeros((243, 243), dtype=np.float64)
        self.ex = np.zeros((243, 243), dtype=np.float64)
        self.ez = np.zeros(972, dtype=np.float64)

        # Variables from other COMMON blocks (ORDER)
        self.fr = np.zeros(140000, dtype=np.float64)
        self.str_val = np.zeros((3, 140000), dtype=np.float64)
        self.bl = np.zeros(140000, dtype=np.float64)
        self.is_val = np.zeros(140000, dtype=np.int32)
        self.j = np.zeros(140000, dtype=np.int32)
        self.iv = np.zeros(140000, dtype=np.int32)

        # Parsed state variables
        self.title = ""
        self.mq = 0
        self.n1 = 0
        self.n2 = 0
        self.nc = 0
        self.iscd = 0
        self.niv = 0
        self.nit = 0
        self.ifpr = 0
        self.nunc = 0
        self.dip = np.zeros(3, dtype=np.float64)
        self.temp = 0.0
        self.kvq = 0

        self.rho1 = 0.0
        self.ivrho1 = 0
        self.beta1 = 0.0
        self.ivbet1 = 0
        self.alpha1 = 0.0
        self.ivalp1 = 0

        self.rho2 = 0.0
        self.ivrho2 = 0
        self.beta2 = 0.0
        self.ivbet2 = 0
        self.alpha2 = 0.0
        self.ivalp2 = 0

        self.jmin = np.zeros(6, dtype=np.int32)
        self.jmax = np.zeros(6, dtype=np.int32)
        self.fmn = np.zeros(6, dtype=np.float64)
        self.fmx = np.zeros(6, dtype=np.float64)
        self.thres = np.zeros(6, dtype=np.float64)

        # Fortran variables A(60,6), IVR(60,6), INPAR(38,6) etc.
        self.a_params = np.zeros((60, 6), dtype=np.float64)
        self.ivr = np.zeros((60, 6), dtype=np.int32)
        self.inpar = np.zeros((38, 6), dtype=np.int32)
        self.ntup = np.zeros(6, dtype=np.int32)
        self.nte = np.zeros(6, dtype=np.int32)

        self.nsig = np.zeros(6, dtype=np.int32)
        self.isig = np.zeros((4, 10, 6), dtype=np.int32)

        # Transitions variables
        self.itra = np.zeros((8, 8191), dtype=np.int32)
        self.ilev = np.zeros((2, 16383), dtype=np.int32)
        self.frq = np.zeros(8191, dtype=np.float64)
        self.wt = np.zeros(8191, dtype=np.float64)
        self.bl_fit = np.zeros(8191, dtype=np.float64)
        self.ntra = 0


def read_input(state, filename):
    """
    Stub for the Fortran INPUT subroutine.
    Reads input parameters and transition data from the given file.
    Updates the shared state.
    """
    with open(filename, 'r') as f:
        # Strip trailing whitespace but keep line structure to handle Fortran's line-based READ(5,*)
        lines = [line.strip() for line in f if line.strip() != '']

    state.title = lines[0]

    def parse_line_floats(line):
        return [float(x.replace('D', 'E').replace('d', 'e')) for x in line.split()]

    line_idx = 1

    def get_parts():
        nonlocal line_idx
        p = lines[line_idx].split()
        line_idx += 1
        return p

    # Global
    parts = get_parts()
    state.iscd, state.nc, state.niv, state.nit, state.ifpr, state.nunc = map(int, parts[:6])
    state.dip[0], state.dip[1], state.dip[2] = map(float, parts[6:9])
    state.temp = float(parts[9])
    state.kvq = 0 # Fortran code does KVQ=0 before the global read

    state.mq = 3 - min(max(abs(state.iscd), 1), 2)

    # Rotor 1
    parts = get_parts()
    state.n1, state.rho1, state.ivrho1 = int(parts[0]), float(parts[1]), int(parts[2])
    state.beta1, state.ivbet1 = float(parts[3]), int(parts[4])
    state.alpha1, state.ivalp1 = float(parts[5]), int(parts[6])

    if state.mq == 2:
        parts = get_parts()
        state.n2, state.rho2, state.ivrho2 = int(parts[0]), float(parts[1]), int(parts[2])
        state.beta2, state.ivbet2 = float(parts[3]), int(parts[4])
        state.alpha2, state.ivalp2 = float(parts[5]), int(parts[6])

    for iv in range(state.niv):
        parts = get_parts()
        state.jmin[iv], state.jmax[iv] = int(parts[0]), int(parts[1])
        state.fmn[iv], state.fmx[iv], state.thres[iv] = map(float, parts[2:5])

        # A array (indices 6 to 20 in Python representing 7 to 21 in Fortran)
        a_arr = []
        while len(a_arr) < 15:
            a_arr.extend(parse_line_floats(lines[line_idx]))
            line_idx += 1
        state.a_params[6:21, iv] = a_arr[:15]

        # IVR array (indices 6 to 20)
        ivr_arr = []
        while len(ivr_arr) < 15:
            ivr_arr.extend([int(x) for x in lines[line_idx].split()])
            line_idx += 1
        state.ivr[6:21, iv] = ivr_arr[:15]

        # Tunneling Params
        param_count = 0
        while True:
            parts = get_parts()
            iq1, iq2, meg = map(int, parts[:3])

            if meg == 0:
                break

            kap, jp, kp = map(int, parts[3:6])
            par = float(parts[6].replace('D', 'E').replace('d', 'e'))
            ivar = int(parts[7])
            scpp = float(parts[8].replace('D', 'E').replace('d', 'e')) if len(parts) > 8 else 1.0

            # Populate a_params and ivr according to Fortran. The logic there:
            # SPAR(KP) = PAR, IGR(KP) = IVAR, SCPP
            # A(KP+21, IV) = SPAR(KP)
            state.a_params[param_count + 21, iv] = par
            state.ivr[param_count + 21, iv] = ivar
            # To mirror Fortran's `INPAR`, we store packed integers. We will just store the individual components for now
            # as that's an implementation detail for Phase 2.

            param_count += 1

        state.ntup[iv] = param_count

        # Symmetry blocks
        sig_count = 0
        while True:
            parts = get_parts()
            is1 = int(parts[0])
            if is1 < 0:
                break
            is2, isp1, isp2 = map(int, parts[1:4])
            state.isig[0, sig_count, iv] = is1
            state.isig[1, sig_count, iv] = is2
            state.isig[2, sig_count, iv] = isp1
            state.isig[3, sig_count, iv] = isp2
            sig_count += 1

        state.nsig[iv] = sig_count

    # Transitions
    lm = 0
    while line_idx < len(lines):
        parts = get_parts()
        is1 = int(parts[0])
        if is1 < 0:
            break

        if state.kvq != 1:
            is2, jq, nq, j_val, nn = map(int, parts[1:6])
            freq = float(parts[6])
            ble, er = map(float, parts[7:9])
            ivq, iv1 = 0, 0
        else:
            is2, ivq, jq, nq, iv1, j_val, nn = map(int, parts[1:8])
            freq, er, ble = map(float, parts[8:11])

        state.itra[0, lm] = is1
        state.itra[1, lm] = is2
        state.frq[lm] = freq
        state.bl_fit[lm] = ble
        # We store the weight as 1/er^2 if er > 0
        if er > 0.0:
            state.wt[lm] = 1.0 / (er * er)
        else:
            state.wt[lm] = 0.0

        lm += 1

    state.ntra = lm

def iterate_fit(state, *args, **kwargs):
    pass

def f_test(state, *args, **kwargs):
    pass

def predict(state, out_f):
    """
    Stub for the Fortran PREDIC subroutine.
    """
    out_f.write("\n PREDICTIONS\n")
    # Physics prediction logic to be added in Phase 2
    pass

def order(state, out_f):
    """
    Stub for the Fortran ORDER subroutine.
    """
    if state.ifpr == 4:
        icatid = 0
        q = 1.0
        logstr0 = -10.0
        logstr1 = -10.0
        out_f.write("\n\n" + "*"*80 + "\n")
        out_f.write(" CATALOG ENTRY INFORMATION\n")
        out_f.write(f" Catalog ID{icatid:21}\n")
        out_f.write(f" PARTITION FUNCTION    {q:22.8e}\n")
        out_f.write(f" 1st INT CUTOFF (LOGSTR0) {logstr0:16.8e}\n")
        out_f.write(f" 2nd INT CUTOFF (LOGSTR1) {logstr1:16.8e}\n")
        out_f.write("*"*80 + "\n\n")

    out_f.write("\n TRANSITIONS ORDERED BY FREQUENCY\n\n")

def main():
    parser = argparse.ArgumentParser(description="ERHAM: Effective rotational Hamiltonian for molecules")
    parser.add_argument("input_file", help="Input file name")
    parser.add_argument("output_file", help="Output file name")
    args = parser.parse_args()

    state = ErhamState()
    now = datetime.now()
    date_str = now.strftime("%Y %b %d %H:%M:%S").lower()
    header = f" program ERHAM Python Port *** date and time: {date_str}\n"

    print("program ERHAM Python Port")
    print(header)

    try:
        with open(args.output_file, 'w') as out_f:
            out_f.write(header)

            read_input(state, args.input_file)

            if state.nit > 0:
                iterate_fit(state)

            f_test(state)
            predict(state, out_f)
            order(state, out_f)

            now = datetime.now()
            date_str = now.strftime("%Y %b %d %H:%M:%S").lower()
            closing = f" program ERHAM Python Port *** date and time: {date_str}\n"
            out_f.write(closing)

    except IOError as e:
        print(f"Error handling files: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
