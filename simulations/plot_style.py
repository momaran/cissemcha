
import os
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="darkgrid")
beta_spins_color = "#b7b1f2"
alpha_spins_color = "#e78895"
beta_spins_color_1 = "#91ddcf"
alpha_spins_color_1 = "#f9b572"

def setup_matplotlib():
    font = {'family' : "Times New Roman", 'weight' : 'bold', 'size'   : 18}
    plt.rc('font', **font)
    plt.rc('legend', fontsize = 10)
    plt.rcParams['lines.linewidth'] = 2.5
    plt.rcParams['lines.markersize'] = 6
    plt.rcParams["figure.figsize"] = [10, 7]
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams['figure.dpi'] = 150

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
