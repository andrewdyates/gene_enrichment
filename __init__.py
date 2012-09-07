#!/usr/bin/python
import random
from enriched import *

import matplotlib 
matplotlib.use('agg') # required for use on OSC servers
import matplotlib.pyplot as plt
import numpy as np

COLORS = ['r', 'b', 'g', 'y', 'c', 'm', 'k', '0.50']
PATTERNS = ['-', ':', '-.']
STYLES = []
for pattern in PATTERNS:
  for color in COLORS:
    STYLES.append({'color': color, 'linestyle': pattern})


def random_indices(n, size):
  return random.sample(xrange(int(n*(n-1)/2)), size)


def make_enrichment_curve_figure(Ranks, title=None, plotpath=None):
  if title is None:
    title = "Gene Enrichment"
  plt.clf(); plt.cla()
  plot_enrichments(Ranks)
  # top left, 8px
  leg = plt.legend(loc=2)
  plt.xlabel('Rank')
  plt.ylabel('# Matched')
  ltext  = leg.get_texts()  # all the text.Text instance in the legend
  plt.setp(ltext, fontsize='xx-small')    # the legend text fontsize
  plt.title(title)
  if plotpath is None:
    plt.show()
  else:
    plt.savefig(plotpath, dpi=600)


def plot_enrichments(Ranks):
  """PyPlot a dictionary of ranks."""
  style_i = 0
  for name, ranks in Ranks.items():
    x = ranks
    y = np.arange(len(x))+1
    plt.plot(x, y, label=name, **STYLES[style_i]);
    style_i += 1
  

                     
