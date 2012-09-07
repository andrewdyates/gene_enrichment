#!/usr/bin/python
"""PINA protein-protein interaction database
http://cbg.garvan.unsw.edu.au/pina/interactome.stat.do


PATH = "/Users/z/Dropbox/archive/Homo_sapiens-20110628.txt"

"ID(s) interactor A"	"ID(s) interactor B"	"Alt. ID(s) interactor A"	"Alt. ID(s) interactor B"	"Alias(es) interactor A"	"Alias(es) interactor B"	"Interaction detection method(s)"	"Publication 1st author(s)"	"Publication Identifier(s)"	"Taxid interactor A"	"Taxid interactor B"	"Interaction type(s)"	"Source database(s)"	"Interaction identifier(s)"	"Confidence value(s)"	"Experimental role(s) interactor A"	"Experimental role(s) interactor B"	"Properties interactor A"	"Properties interactor B"	"HostOrganism(s)"

Example interaction column entry
'MI:0007(anti tag coimmunoprecipitation)|MI:0007(anti tag coimmunoprecipitation)'
"""
import re
from enriched import *
from lab_util.gene_clean import *
import sys

RX_INTERACT = re.compile('[^(]+([^)]+)\)')

class PINAEnriched(Enriched):
  """Enrichment object for PINA protein interaction list.
  Cleans variable names before comparing.
  """
  RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")
  NAME_COL_A = 2
  NAME_COL_B = 3

  def __init__(self, fp):
    """Initialize from new open fp to PINA formatted file."""
    self.headers = fp.next().rstrip('\n').split('\t')
    self.pair_type = {}
    super(PINAEnriched, self).__init__(pairs=self._pair_gen(fp))
  
  def _pair_gen(self, fp):
    for line in fp:
      row = line[:-1].split('\t')
      s_a, s_b = (row[self.NAME_COL_A], row[self.NAME_COL_B])
      m_a, m_b = self.RX_GENE_NAME.match(s_a), self.RX_GENE_NAME.match(s_b)
      if not (m_b and m_b):
        continue
      pair = (clean(m_a.group(1)), clean(m_b.group(1)))
      self._add_pair_type(pair, row)
      yield pair

  def _hash(self, x,y):
    return ",".join(sorted((x,y)))
  
  def _add_pair_type(self, pair, row):
    hash = self._hash(*pair)
    # Corresponds to "Interaction detection method(s)"
    s = row[6]
    m = RX_INTERACT.match(s)
    if not m:
      print "!", s
      sys.exit(1)
    self.pair_type[hash] = m.group(1)
    
  def get_relation_types(self, x, y):
    """Return list of relationship types for gene pair x, y."""
    if not self.exists(x,y):
      return None

    

  def exists(self, x, y):
    return super(PINAEnriched, self).exists(clean(x), clean(y))
  def is_in(self, x):
    return super(PINAEnriched, self).is_in(clean(x))

  def pair_indices(self, varlist):
    cleaned_varlist = [clean(s) for s in varlist]
    return super(PINAEnriched, self).pair_indices(cleaned_varlist)

  def var_indices(self, varlist):
    cleaned_varlist = [clean(s) for s in varlist]
    return super(PINAEnriched, self).var_indices(cleaned_varlist)
  
