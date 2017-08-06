import pandas as pd
results = pd.read_csv("/Users/darrin/git/cutlery/piNpiS_random_simulator/test_results.csv")

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter
import numpy as np
import sys
import os

Rf = 65 / 255
Bf = 85 / 255
pdict = {'red': ((0.0, Rf, Rf),
                 (1.0, Rf, Rf)),
         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, Bf, Bf),
                  (1.0, Bf, Bf)),
         'alpha': ((0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0))
         }
# Now we will use this example to illustrate 3 ways of
# handling custom colormaps.
# First, the most direct and explicit:
purple1 = LinearSegmentedColormap('Purple1', pdict)
    
def generate_heat_map(panel, data_frame, color):
    # This single line controls plotting the hex bins in the panel
    hex_vals = panel.hexbin(data_frame['pi'], data_frame['piNpiS'], gridsize=10,
                            linewidths=0.0, cmap=color)
    #for each in panel.spines:
    #    panel.spines[each].set_visible(False)

    counts = hex_vals.get_array()
    return counts

plt.style.use('BME163')

#set the figure dimensions
figWidth = 5
figHeight = 5
plt.figure(figsize=(figWidth,figHeight))

#set the panel dimensions
panelWidth = 3
panelHeight = 3

#find the margins to center the panel in figure
leftMargin = (figWidth - panelWidth)/2
bottomMargin = ((figHeight - panelHeight)/2) + 0.25

panel0=plt.axes([leftMargin/figWidth, #left
                 bottomMargin/figHeight,    #bottom
                 panelWidth/figWidth,   #width
                 panelHeight/figHeight])     #height
panel0.tick_params(axis='both',which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on', \
                   right='off', labelright='off',\
                   top='off', labeltop='off')               

panel0.set_xlim([min(results['pi'])*0.9, max(results['pi'])*1.1])
panel0.set_ylim([0, max(results['piNpiS']) * 1.1])
panel0.plot(results['pi'], results['piNpiS'], 'bo', alpha = 0.1, ms = 2)
#counts = generate_heat_map(panel0, results, purple1)
plt.show()