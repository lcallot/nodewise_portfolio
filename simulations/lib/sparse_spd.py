# python script to generate sparse spd matrices
# called by rpython

import numpy as np
from sklearn.datasets.samples_generator import make_sparse_spd_matrix


def mk_spd(dim,alpha,maxc,minc,rs):
  prec = make_sparse_spd_matrix(dim=dim, alpha=alpha, largest_coef=maxc, smallest_coef=minc, random_state=rs)
  prec = prec.tolist()
  return[prec]
