from __future__ import print_function
import rebound
from scipy.stats import kstest, ks_2samp
import numpy as np
from astropy.io import ascii
import sys

@np.vectorize
def ethermal(ecc):
	if ecc<0:
		return 0
	else:
		return min(ecc**2.0,1.0)

loc='./'
base=sys.argv[1]
ff='trial_'+base+'.csv'
dat=ascii.read(loc+ff, fill_values=(('I', np.nan), ('-I', np.nan), ('Coll', np.nan)))
ff='trial_end_'+base+'.csv'
dat_final_pos=ascii.read(loc+ff, names=['t', 'x', 'y','vx', 'vy'], data_start=1, delimiter=',')
# print(dat_final_pos[dat_final_pos['t']<0.999])

sim=rebound.Simulation()
sim.G=1


filt=~dat['as/abin'].mask
f=open('trial_end_{0}_orb_analysis.csv'.format(base), 'w')
f.write('row, a, e, ai, ei, S\n')

for idx in range(len(dat)):
	##Consider only disrupted binaries, and discard cases in which the integration stopped short.
	if filt[idx] or dat_final_pos['t'][idx]<0.999:
		continue
	qq=dat[idx]['q']
	mm=dat[idx]['m']
	sim.add(m=mm/(1.0+qq), x=0, y=0, z=0, vx=0, vy=0, vz=0)
	sim.add(m=mm*qq/(1.0+qq),x=dat_final_pos[idx]['x'], y=dat_final_pos[idx]['y'], z=0,\
		vx=dat_final_pos[idx]['vx'], vy=dat_final_pos[idx]['vy'], vz=0)

	orb=sim.particles[0].calculate_orbit(primary=sim.particles[1])
	f.write(str([idx, orb.a, orb.e, dat[idx]['abin'],\
				 dat[idx]['e'], dat[idx]['S']]).replace('[','').replace(']',''))
	f.write('\n')
	sim.remove(1)
	sim.remove(0)
f.close()
