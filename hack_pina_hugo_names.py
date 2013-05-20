"""Read PINA file, print simple adjaceny list"""
from pina import *

pina_fpath = "Homo sapiens-20121210.txt"
FNAME_HGNC = "/Users/z/Dropbox/biostat/git_repos/transcription_factors/hgnc_alias_list.txt"
P = PINAEnriched(open(pina_fpath))
row_names, A = P.return_adj_matrix()

syms = {}
official= set()
fp = open(FNAME_HGNC)
fp.next()
bad = 0
for line in fp:
  row = line.strip().split('\t')
  if len(row) < 6:
    continue
  sym, aliases = row[1], filter(None, row[4].split(', ')+row[5].split(', '))
  syms[sym] = sym
  official.add(sym)
  for s in aliases:
    if s in syms:
      bad += 1
    syms[s] = sym # will override dupes... oh well
print len(syms)
print bad


for s in row_names:
  if s not in official:
    print "Not official:", s
    if s in syms:
      print "  Replacement:", syms[s]
