"""
Friends-of-Friends (FOF) for N-body simulations

Peter Creasey - Oct 2016

"""
from __future__ import absolute_import, print_function
from hfof.lib import fof3d, get_cells
from lizard.periodic import pad_unitcube
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix, csgraph
from numpy import square, flatnonzero, ones, zeros_like, cumsum, concatenate, \
    arange, searchsorted, bincount, sort, diff, int8, argsort, array, unique, empty, int32, uint32


def bin_id_data(ids, log):
    """
    Helper-function to bin up the ids. 
    ids - the ID for N items
    log - write out debug data

    returns sort_idx, bin_data

    s.t.
    sort_idx - an array such that ids[sort_idx] is sorted
    bin_data - a (Q,3) array such that Q is the number of unique IDs, 
        bin_data[i,0] is the count of the i-th smallest ID and
        ids[sort_idx][bin_data[i,1]:bin_data[i,2]] would all be the i-th
        smallest ID.
               
    Since these indices are used in the C-function, in principle you can cause 
    a seg-fault if you get them wrong. Mess with this function at your own peril!
    """

    num_pos = len(ids)
    print('Finding unique values', file=log)
    filled_cells, bin_fill = unique(ids, return_counts=True)

    num_bins = len(filled_cells)
    av_bin_fill = float(num_pos)/num_bins
    av_fill_per_id = float(square(bin_fill).sum())/num_pos
    print('Put {:,} IDs into {:,} bins'.format(num_pos, num_bins),
          'of fill {:,}-{:,} IDs,'.format(bin_fill.min(), bin_fill.max()), 
          'average %.2f.'%av_bin_fill,
          'Average ID lives in a bin of fill %.2f IDs'%av_fill_per_id, file=log)


    


    bin_data = empty((num_bins,2), dtype=int32)
    bin_data[:,1]   = cumsum(bin_fill)
    bin_data[:,0]   = bin_data[:,1] - bin_fill
    del bin_fill
    print('Sorting {:,} IDS'.format(num_pos), file=log)
    sort_idx = argsort(ids).astype(uint32) 
    del ids
    return sort_idx,(filled_cells, bin_data)

