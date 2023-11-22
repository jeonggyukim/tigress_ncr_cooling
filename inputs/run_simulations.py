#!/usr/bin/env python

import numpy as np
import re
import subprocess
import argparse
import shutil
import os

parser = argparse.ArgumentParser()
parser.add_argument('--shld', action='store_true', default=False,
                    help='Turn on shielding')
parser.add_argument('--rundir', type=str, default='../bin',
                    help='Directory where simulations are run')
args = vars(parser.parse_args())

locals().update(args)

if __name__ == '__main__':

    if not os.path.exists(rundir):
        os.makedirs(rundir)

    nomake = True
    tncinput_template = 'template.tncinput'

    t_end_Myr = 1e3

    xi_cr0_def = 2e-16

    # Log gas metallicity, dust abundance, and radiation
    # z_g = np.arange(-3.0,1.0,0.5)
    z_g = np.array([-2.0,-1.0,0.0, np.log10(3.0)])
    z_d = z_g

    suite = 'fiducial'
    if suite == 'fiducial':
        # Case A. Make radiation and CR scale together
        chi0 = np.arange(-2.0,3.0,1.0)
        xi_cr0 = np.arange(-2.0,3.0,1.0)
    elif suite == 'var_cr':
        # Case B. Fix radiation field but change only CRs
        xi_cr0 = np.arange(-2.0,1.0,0.25)
        chi0 = np.repeat(0.0,len(xi_cr0))

    if shld:
        len_shld0 = 5.0
        flag_cool_dust = 1
        shld_str = 'shld'
    else:
        len_shld0 = 0.0
        flag_cool_dust = 0
        shld_str = 'noshld'

    z = (z_g,z_d)
    chi = (chi0,xi_cr0)
    for z_g_,z_d_ in zip(z_g,z_d):
        for chi0_,xi_cr0_ in zip(chi0,xi_cr0):
            print('z_g,z_d,chi0,xi_cr0',z_g_,z_d_,chi0_,xi_cr0_)
            if suite == 'fiducial':
                problem_id = '{0:s}_zg{1:.2f}_zd{2:.2f}_chi{3:.1f}'.format(
                    shld_str, z_g_, z_d_, chi0_)
            elif suite == 'var_cr':
                problem_id = '{0:s}_zg{1:.1f}_zd{2:.1f}_xi{3:.2f}'.format(
                    shld_str, z_g_, z_d_, xi_cr0_)

            tncinput_out = 'tncinput.{:s}'.format(problem_id)

            subs = dict()
            subs['problem_id'] = problem_id
            subs['z_g'] = str(10.0**z_g_)
            subs['z_d'] = str(10.0**z_d_)
            subs['chi0'] = str(10.0**chi0_)
            subs['xi_cr0'] = str(xi_cr0_def*10.0**xi_cr0_)

            subs['len_shld0'] = str(len_shld0)
            subs['flag_cool_dust'] = str(flag_cool_dust)
            subs['t_end_Myr'] = str(t_end_Myr)

            # Read template, make substitutions, and write a file
            with open(tncinput_template, 'r') as f:
                tmp = f.read()
            for k, v in subs.items():
                tmp = re.sub(r'@{0}@'.format(k), v, tmp)
            with open(tncinput_out, 'w') as f:
                f.write(tmp)

            # Run simulations
            shutil.move(tncinput_out, os.path.join(rundir, tncinput_out))
            cmd = ['../bin/tigress_ncr_cooling',
                   '-i', os.path.join(rundir, tncinput_out),
                   '-d', rundir]

            out = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
            print(out.decode('unicode_escape'))
            print('')
