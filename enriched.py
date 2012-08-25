#!/usr/bin/python
from py_symmetric_matrix import *

class Enriched(object):
  """Enrichment pair set. Intented for use with gene enrichment, but could be used for any ranking."""

  def __init__(self, pairs):
    """Initialize.
    Args:
      pairs [*(str, str)] iterable of pairs of names.
    """
    self.vars = set()
    self.n_pairs = 0
    # if A < B, then A is key and B is in set `self.d_pairs[A]`
    self.d_pairs = {}
    
    for x, y in pairs:
      x, y = sorted((x,y))
      if not (x and y) or x == y:
        continue
      if not (x in self.d_pairs and y in self.d_pairs[x]):
        self.n_pairs += 1
        self.d_pairs.setdefault(x, set()).add(y)
        self.vars.add(x)
      
  def exists(self, x, y):
    """Return if (x, y) variables are a pair in this Enrichment set."""
    x, y = sorted((x,y))
    if not (x and y):
      return False
    if not (x in self.vars and y in self.vars):
      return False
    return y in self.d_pairs[x]

  def pair_indices(self, varlist):
    """Return indices of self.pairs for sym dependency matrix indexed by varlist.

    Args:
      varlist: [str] of ordered variable names
    Returns:
      [int] of indices to variable pairs in squareform all-pairs matrix with `varlist` labeled rows.
    """
    n = len(varlist)
    idxs = []
    var_d = dict([(v, i) for i, v in enumerate(varlist)])
    for x, adj_set in self.d_pairs.items():
      for y in adj_set:
        try:
          i, j = var_d[x], var_d[y]
        except KeyError:
          continue
        idxs.append(sym_idx(i,j,n))
    return idxs

  def var_indices(self, varlist):
    """Return row indices of varlist for symbols in enrichment set.

    Args:
      varlist: [str] of ordered variable names
    Returns:
      [int] of indices in varlist & self.geneset
    """
    idxs = []
    for i, v in enumerate(varlist):
      if v in self.vars:
        idxs.append(i)
    return idxs
