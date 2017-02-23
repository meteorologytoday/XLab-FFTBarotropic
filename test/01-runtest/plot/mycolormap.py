import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


gray = [(0.0, 1.0, 1.0),
        (0.15, 0.5, 0.5),
        (1.0,0.0,0.0)]

cdict = {'red': gray, 'green': gray, 'blue': gray}

cmap_vorticity = LinearSegmentedColormap('vorticity', cdict)
