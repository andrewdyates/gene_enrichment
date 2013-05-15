"""Read PINA file, print simple adjaceny list"""
from pina import *
import matrix_io as mio

pina_fpath = "/nfs/01/osu6683/Homo sapiens-20121210.txt"
P = PINAEnriched(open(pina_fpath))
row_names, A = P.return_adj_matrix()
mio.save(M=A, fp="/nfs/01/osu6683/pina_adj_matrix.tab", fmt="%d", row_ids=row_names, col_ids=row_names, fill_upper_left=False)
