import numpy as np
import sys, getopt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.pyplot as pplt
import matplotlib.font_manager as fm
import matplotlib as mplt
from matplotlib.patches import Rectangle
from mycolormap import cmap_vorticity


rc('font', **{'family':'sans-serif', 'serif': 'Bitstream Vera Serif', 'sans-serif': 'MS Reference Sans Serif', 'size': 20.0, 'weight' : 100})
rc('axes', **{'labelsize': 20.0, 'labelweight': 100})
rc('mathtext', **{'fontset':'stixsans'})

default_linewidth = 2.0
default_ticksize = 10.0

mplt.rcParams['lines.linewidth'] =   default_linewidth
mplt.rcParams['axes.linewidth'] =    default_linewidth
mplt.rcParams['xtick.major.size'] =  default_ticksize
mplt.rcParams['xtick.major.width'] = default_linewidth
mplt.rcParams['ytick.major.size'] =  default_ticksize
mplt.rcParams['ytick.major.width'] = default_linewidth
mplt.rcParams['xtick.major.pad']='12'
mplt.rcParams['ytick.major.pad']='12'


import config as cf
from mycolormap import cmap_vorticity


try:
	opts, args = getopt.getopt(sys.argv[1:], "", ["input-dir=", "output-dir=", "start-step="])
except getopt.GetoptError as err:
	print(err)
	sys.exit(2)

in_dir = None
out_dir = None
start_step = 0
barb_skip = cf.barb_skip
for o, a in opts:
	if o == "--input-dir":
		in_dir = a
	elif o == "--output-dir":
		out_dir = a
	elif o == "--start-step":
		start_step = int(a)

if in_dir is None or out_dir is None:
	print("Must give paramter --input-dir and --output-dir")
	sys.exit(2)

print("Input directory  : %s" % (in_dir,))
print("Output directory : %s" % (out_dir,))
print("Start step       : %d" % (start_step,))

x_vec = cf.x_vec / 1000.0
y_vec = cf.y_vec / 1000.0


cb_vec = np.linspace(0, 20, num = 21)

pres_interval = np.arange(-100, 10, 5)

# calculate figure size
graph_size = np.array([4.0, 3.0]) * 3.0;
space = {'wspace': 2.0, 'hspace': 1.8, 'tspace': 1.0, 'bspace': 1.0, 'lspace': 1.0, 'rspace': 2.5, 'break':0.1}

figsize = [graph_size[0] + space['lspace'] + space['rspace'],
           graph_size[1] + space['tspace'] + space['bspace']];

figdpi = 100



fig_info = {
	'left'  : space['lspace'] / figsize[0],
	'right' : 1.0 - (space['rspace'] / figsize[0]),
	'bottom': space['bspace'] / figsize[1],
	'top'   : 1.0 - (space['tspace'] / figsize[1]),
	'wspace': space['wspace'] / graph_size[0],
	'hspace': space['hspace'] / graph_size[1]
};

graph_rect = (
		space['lspace'] / figsize[0],
		space['bspace'] / figsize[1],
		graph_size[0] / figsize[0],
		graph_size[1] / figsize[1]
);

t_cnt = 0.0

for step in range(int(start_step / cf.record_step) * cf.record_step, cf.total_steps, cf.record_step):

	i = int(step / cf.record_step)
	try:
		vort = np.fromfile("%s/vort_step_%d.bin" % (in_dir, step), dtype='<f4', count=(cf.nx*cf.ny)).reshape((cf.nx, cf.ny)).transpose()
		u = np.fromfile("%s/u_step_%d.bin" % (in_dir, step), dtype='<f4', count=(cf.nx*cf.ny)).reshape((cf.nx, cf.ny)).transpose()
		v = np.fromfile("%s/v_step_%d.bin" % (in_dir, step), dtype='<f4', count=(cf.nx*cf.ny)).reshape((cf.nx, cf.ny)).transpose()
	
	except IOError as e:
		print(e)
		print(e.strerror)
		raise Exception("Error when loading files... quit program.")
	

	fig = pplt.figure(figsize=figsize);
	ax = fig.add_axes(graph_rect, autoscale_on=False);
	ax.set_xlim([x_vec[0], x_vec[-1]])
	ax.set_ylim([y_vec[0], y_vec[-1]])
	ax.set_xlabel(r'x [$\mathrm{km}$]')
	ax.set_ylabel(r'y [$\mathrm{km}$]')
	ax.set_aspect(1)

	cax = fig.add_axes([0.85, 0.1, 0.05,0.8])



	cbar_mappable = ax.contourf(x_vec, y_vec, vort * 1000.0, cb_vec, cmap=cmap_vorticity)#plt.get_cmap("gray_r"))
	cbar = fig.colorbar(cbar_mappable, cax=cax, orientation='vertical')
	#strm = ax.streamplot(x_vec, y_vec, u, v, linewidth=2)

	ax.barbs(x_vec[::barb_skip], y_vec[::barb_skip], u[::barb_skip,::barb_skip] * 0.5144, v[::barb_skip,::barb_skip] * 0.5144, length=8)
	ax.text(1.1, 0.5, r'$\zeta$ [$\times\,10^{-3}\,\mathrm{s}^{-1}$]', rotation=90, horizontalalignment='left', verticalalignment="center", transform=ax.transAxes, fontsize=30)
	

	ax.text(10, 10, "%02d:%02d:%02d" % (int(t_cnt/3600), int(t_cnt/60) % 60, t_cnt % 60))
	t_cnt = step * cf.dt

	file_out = "%s/step_%d.png" % (out_dir, step) 
	fig.savefig(file_out, dpi=figdpi, format='png')
	print("Output image: %s" % (file_out,))
	pplt.close(fig)
