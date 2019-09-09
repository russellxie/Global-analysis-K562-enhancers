#!/usr/bin/env python3
import collections
import tables
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from multiprocessing import Pool


def _find_non_zero_cells(df, gbm):
    cell_bc = df.T.loc[df.sum() != 0].index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(gbm.barcodes == i.encode('utf-8'))
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]
