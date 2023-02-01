#!/usr/bin/env python

import numpy as np
import re
import subprocess
import argparse
import shutil
import os.path as osp

parser = argparse.ArgumentParser()
parser.add_argument('--shld', type=bool, default=False,
                    help='Turn on shielding')
parser.add_argument('--rundir', type=str, default='../bin',
                    help='Directory where simulations are run')
args = vars(parser.parse_args())

locals().update(args)

if __name__ == '__main__':

    nomake = True
    tncinput_template = 'template.tncinput'
    
    t_end_Myr = 1e3

    # Log gas metallicity and dust abundance
    z_g = np.arange(-3.0,1.0,0.5)
    z_d = z_g
    chi0 = 1.0
    xi_cr0 = 2e-16

    if shld:
        len_shld0 = 5.0
        flag_dust_cool = 1
        shld_str = 'shld'
    else:
        len_shld0 = 0.0
        flag_dust_cool = 0
        shld_str = 'noshld'

    for i,(z_g_,z_d_) in enumerate(zip(z_g,z_d)):
        print(z_g_,z_d_)
        problem_id = '{0:s}_zg{1:.1f}_zd{2:.1f}_chi{3:.1f}'.format(
            shld_str, z_g_, z_d_, chi0)
        tncinput_out = 'tncinput.{:s}'.format(problem_id)
        
        subs = dict()
        subs['problem_id'] = problem_id
        subs['z_g'] = str(10.0**z_g_)
        subs['z_d'] = str(10.0**z_d_)
        subs['xi_cr0'] = str(xi_cr0)
        subs['chi0'] = str(chi0)

        subs['len_shld0'] = str(len_shld0)
        subs['flag_dust_cool'] = str(flag_dust_cool)
        subs['t_end_Myr'] = str(t_end_Myr)
        
        # Read template, make substitutions, and write a file
        with open(tncinput_template, 'r') as f:
            tmp = f.read()
        for k, v in subs.items():
            tmp = re.sub(r'@{0}@'.format(k), v, tmp)
        with open(tncinput_out, 'w') as f:
            f.write(tmp)

        # Run simulations
        shutil.move(tncinput_out, osp.join(rundir, tncinput_out))
        cmd = ['../bin/tigress_ncr_cooling',
               '-i', osp.join(rundir, tncinput_out),
               '-d', rundir]
        
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
        print(out.decode('unicode_escape'))
        print('')

