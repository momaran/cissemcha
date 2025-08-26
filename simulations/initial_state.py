#--------------------------------------
#Libreries to be imported:
import matplotlib.pyplot as plt
from matplotlib import rc
from PIL import Image
import imageio
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.animation as animation
#--------------------------------------

# Setting Plot style:
sns.set_theme(style="darkgrid")

# Option 0:
beta_spins_color = "#b7b1f2"
alpha_spins_color = "#e78895"

# Option 1:
beta_spins_color_1 = "#91ddcf"
alpha_spins_color_1 = "#f9b572"

#-------------------------------------------------
#       Setting configurations
#-------------------------------------------------

font = {'family' : "Times New Roman",
       'weight' : 'bold',
       'size'   : 15  }

plt.rc('font', **font)
plt.rc('legend', fontsize = 10) # using a size in points
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['lines.markersize'] = 6
plt.rcParams["figure.figsize"] = [10, 6]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['figure.dpi'] = 150  # Increase display resolution

