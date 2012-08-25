#!/usr/bin/python
import random
from enriched import *

def random_indices(n, size):
  return random.sample(xrange(int(n*(n-1)/2)), size)


