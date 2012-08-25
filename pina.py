#!/usr/bin/python
"""PINA protein-protein interaction database
http://cbg.garvan.unsw.edu.au/pina/interactome.stat.do
"""
import re
from enriched import *
from lab_util.gene_clean import *

class PINAEnriched(Enriched):
  """Enrichment object for PINA protein interaction list.
  Cleans variable names before comparing.
  """
  RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")
  NAME_COL_A = 2
  NAME_COL_B = 3

  def __init__(self, fp):
    """Initialize from new open fp to PINA formatted file."""
    super(PINAEnriched, self).__init__(pairs=self._pair_gen(fp))
  
  def _pair_gen(self, fp):
    for line in fp:
      row = line[:-1].split('\t')
      s_a, s_b = (row[self.NAME_COL_A], row[self.NAME_COL_B])
      m_a, m_b = self.RX_GENE_NAME.match(s_a), self.RX_GENE_NAME.match(s_b)
      if not (m_b and m_b):
        continue
      yield (clean(m_a.group(1)), clean(m_b.group(1)))

  def exists(self, x, y):
    return super(PINAEnriched, self).exists(clean(x), clean(y))

  def pair_indices(self, varlist):
    cleaned_varlist = [clean(s) for s in varlist]
    return super(PINAEnriched, self).pair_indices(cleaned_varlist)

  def var_indices(self, varlist):
    cleaned_varlist = [clean(s) for s in varlist]
    return super(PINAEnriched, self).var_indices(cleaned_varlist)
  

