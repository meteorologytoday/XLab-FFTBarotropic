import numpy as np;

nx = 768
ny = 768

Lx = 600000.0
Ly = 600000.0

x_range = np.array([0.0, Lx])
y_range = np.array([0.0, Ly])

x_vec = np.linspace(0, Lx, nx)
y_vec = np.linspace(0, Ly, ny)

t_range = np.array([0, 12.0])

dt = 3.0
record_step = 100
total_steps = int(12 * 60 * 60 / dt)
avg_n = 30

barb_skip=48
