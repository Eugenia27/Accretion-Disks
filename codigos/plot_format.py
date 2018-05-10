import matplotlib as mpl 
import matplotlib.pyplot as plt

def plot_format(): 
    #plt.rc('text', usetex=True)
    font = {'family': 'sans-serif', 'size': 33, 'serif': ['computer modern roman']}
    plt.rc('font', **font)
    plt.rc('legend', **{'fontsize': 24}) 
    #plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['xtick.major.size'] = 8
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['ytick.major.size'] = 6
    plt.rcParams['ytick.minor.size'] = 3
    plt.rc('lines', linewidth=3)


