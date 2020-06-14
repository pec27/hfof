"""
Tools for making the neighbour values (masks and hashed walk values)
"""
from __future__ import print_function
from numpy import mgrid, argsort, cumsum

def ngb_octants():
    z,y,x = mgrid[:3,:3,:3]-1
    x,y,z = x.ravel(), y.ravel(), z.ravel()

    order = argsort(x*x+y*y+z*z)
    walk_ngbs = []
    hash_walk = []
    octant_masks = [[] for i in range(8)]
    for i,j,k in zip(x[order],y[order], z[order]):
        if i<0:
            continue
        elif i==0:
            if j<0:
                continue
            elif j==0:
                if k<=0:
                    continue

        s = ''
        for cpt in (['','M'][i], ['-N','','N'][j+1], ['-1','','1'][k+1]):
            if len(s) and len(cpt)==1:
                s = '+'.join((s,cpt))
            else:
                s = s + cpt
        walk_ngbs.append(s)
        if s=='1':
            hash_walk.append('hprime&hmask')
        elif len(s)==1:
            hash_walk.append('%s*hprime&hmask'%s)
        else:
            hash_walk.append('(%s)*hprime&hmask'%s)

        # Loop over octants
        for octant in range(8):
            oi, oj, ok = octant>>2, (octant>>1)&1, octant&1
            mask = 0 # 8 bits
            for adj_octant in range(8):
                ai, aj, ak = adj_octant>>2, (adj_octant>>1)&1, adj_octant&1
                # relative separation of octants
                sep = 2*i+oi-ai, 2*j+oj-aj, 2*k+ok-ak
                mdist2 = sum(max(abs(v)-1,0) for v in sep)
                if mdist2<3:
                    # octant width is b/sqrt(3)
                    mask |= 1<<adj_octant

            octant_masks[octant].append(mask)
#            print octant, mask, hex(mask)
            
        print(i,j,k, s)
    print('const unsigned int walk_ngbs[13] = {'+', '.join(walk_ngbs)+'};')
    print('const unsigned int hash_ngb[13] = {'+', '.join(hash_walk)+'};')


    print('const unsigned char octant_masks[104] = {'+',\n'.join(', '.join(hex(x) for x in o) for o in octant_masks)+'};')


def disp_corner(sc, L, ct=4):
    """
    Nearest displacement along blocks of 4
    """
    return min(abs((x-L*ct) - sc) for x in range(ct))

def _unpack_idx(idx, nbits, dim):
    """ Convert from raster coords to (i_1, i_2,..., i_dim) """
    mask = (2**nbits)-1
    res = []
    for i in range(dim):
        res.append((idx>>(nbits*(dim - 1 - i))) & mask)
    return tuple(res)

def find_andmask(ijk, sc, nbits=2, dim=3):
    """ Find which cells in block ijk are visible from sc """
    o_ijk = _unpack_idx(sc, nbits, dim)

    adj_idx = []
    all_mask = 0
    tot_bits = dim * nbits 
    tot_cells = 2**tot_bits

    block_width = 2**nbits

    for adj in range(tot_cells):
        adj_ijk = _unpack_idx(adj, nbits, dim)

        sep = [abs(o - (a-block_width*p)) for o,a,p in zip(o_ijk, adj_ijk, ijk)]

        min_dist2 = sum(max(x-1,0)**2 for x in sep)
        # Is this cell visible from o_ijk?
        if min_dist2<dim:
            adj_idx.append(adj)
            all_mask |= 1<<adj

    # 6 different bits to check
    mask = 0
    for i in range(tot_bits):
        a1 = (adj_idx[0]>>i)&1
        if all((adj>>i)&1==a1 for adj in adj_idx):
            mask |= 1<<i
    
    return (mask, adj_idx[0]&mask, all_mask)


def ngb_64():
    z,y,x = mgrid[:3,:3,:3]-1
    x,y,z = x.ravel(), y.ravel(), z.ravel()

    order = argsort(x*x+y*y+z*z)
    walk_ngbs = []
    hash_walk = []

    adj_walk = []
    for i,j,k in zip(x[order],y[order], z[order]):
        if i<0:
            continue
        elif i==0:
            if j<0:
                continue
            elif j==0:
                if k<=0:
                    continue

        s = ''
        for cpt in (['','M'][i], ['-N','','N'][j+1], ['-1','','1'][k+1]):
            if len(s) and len(cpt)==1:
                s = '+'.join((s,cpt))
            else:
                s = s + cpt
        walk_ngbs.append(s)
        if s=='1':
            hash_walk.append('hprime&hmask')
        elif len(s)==1:
            hash_walk.append('%s*hprime&hmask'%s)
        else:
            hash_walk.append('(%s)*hprime&hmask'%s)

        adj_walk.append((i,j,k))
    quad_lists = []
    all_checks = []
    all_masks=[]
    bit_masks=[]
    # Loop over sub-cells, determine which are possible
    for sc in range(64):
        oi, oj, ok = (sc>>4)&3, (sc>>2)&3, sc&3
        my_ngb = []
        mdists = []
        and_masks = []

        for e,(i,j,k) in enumerate(adj_walk):
            # relative separation of sub-cell

            sep = [disp_corner(oi,i), disp_corner(oj,j), disp_corner(ok,k)]
#            print i,j,k, sep
            mdist2 = sum(max(v-1,0)**2 for v in sep)
            if mdist2<3:
                # octant width is b/sqrt(3)
                my_ngb.append(e)
                mdists.append(mdist2)
                and_masks.append(find_andmask((i,j,k), sc))
        # make magic numbers
        mask = 0
        for i,ngb in enumerate([my_ngb[j] for j in argsort(mdists)]):
            mask |= (ngb&0xF)<<(4*i)
        mask |= 13 << (4*len(my_ngb))

        res = []

        bit_masks.append(find_andmask((0,0,0), sc)[2])
        bit_masks.append(mask) # insert the walk
        for i,(ngb,amask,mres,bmask) in enumerate([(my_ngb[j], and_masks[j][0], and_masks[j][1], and_masks[j][2]) for j in argsort(mdists)]):
            res.append(ngb | (amask<<4) | (mres<<10))
            bit_masks.append(bmask)

        # end point 

        quad_lists.append(mask)
        all_checks.append(len(my_ngb)+2)
        print(sc, ', '.join(hex(x) for x in res))
        all_masks.extend(res)

    print('const unsigned int walk_ngbs[13] = {'+', '.join(walk_ngbs)+'};')
    print('const unsigned int hash_ngb[13] = {'+', '.join(hash_walk)+'};')
    hw2 = [i for sub in zip(hash_walk, walk_ngbs) for i in sub]
    print('const unsigned int hash_walk[26] = {'+', '.join(hw2)+'};')
    print('const unsigned int quad_masks[64] = {\n'+', '.join(hex(q) for q in quad_lists)+'};')
#    print('const unsigned char octant_masks[104] = {'+',\n'.join(', '.join(hex(x) for x in o) for o in octant_masks)+'};'
    print('Total walks', sum(all_checks))
    print('static const unsigned char ngb_start[65] = {\n'+', '.join('%d'%v for v in [0] + list(cumsum(all_checks))) + '};\n')
    print('static const uint16_t ngb_mask[220] = {\n'+', '.join(hex(m) for m in all_masks) + '};\n')
    print('static const uint64_t ngb_bits[%d] = {\n'%len(bit_masks)+', '.join(hex(m) for m in bit_masks) + '};\n')
    print('static const uint64_t self_bits[64] = {\n'+', '.join(hex(m) for m in [find_andmask((0,0,0), j)[2] for j in range(64)]) + '};\n')



def ngb2d_64():
    """ 
    like ngb_64, but for the 8x8 cells in 2d 
    Note there are only 4 possible neighbour blocks, and perhaps only 5/6(?) walks?
    """
    nbits, dim = 3,2
    # Searches of adjacent blocks
    y,x = mgrid[:3,:3]-1
    x,y = x.ravel(), y.ravel()

    order = argsort(x*x+y*y)
    walk_ngbs = []
    hash_walk = []

    adj_walk = []
    for i,j in zip(x[order],y[order]):
        if i<0:
            continue
        elif i==0:
            if j<=0:
                continue

        s = ''
        for cpt in (['','N'][i], ['-1','','1'][j+1]):
            if len(s) and len(cpt)==1:
                s = '+'.join((s,cpt))
            else:
                s = s + cpt
        walk_ngbs.append(s)
        if s=='1':
            hash_walk.append('hprime&hmask')
        elif len(s)==1:
            hash_walk.append('%s*hprime&hmask'%s)
        else:
            hash_walk.append('(%s)*hprime&hmask'%s)

        adj_walk.append((i,j))
    quad_lists = []
    all_checks = []
    all_masks=[]
    bit_masks=[]
    # Loop over sub-cells, determine which are possible

    tot_cells = 2**(nbits * dim)
    self_block = (0,0)
    
    assert(tot_cells==64)
    for sc in range(tot_cells):
        o_ij = _unpack_idx(sc, nbits, dim)
        my_ngb = []
        mdists = []
        and_masks = []

        for e,(i,j) in enumerate(adj_walk):
            # relative separation of sub-cell

            sep = [disp_corner(oi,block_i,8) for oi,block_i in zip(o_ij, (i,j))]
#            print('Separation', i,j, sep)
            mdist2 = sum(max(v-1,0)**2 for v in sep)
            if mdist2<2:
                # octant width is b/sqrt(dim)
                my_ngb.append(e)
                mdists.append(mdist2)
                and_masks.append(find_andmask((i,j), sc, nbits, dim))
        # make magic numbers
        mask = 0
        for i,ngb in enumerate([my_ngb[j] for j in argsort(mdists)]):
            mask |= (ngb&0xF)<<(4*i)
        mask |= 13 << (4*len(my_ngb))

        res = []
        # self mask

        #bit_masks.append(find_andmask(self_block, sc, nbits, dim)[2])
        #bit_masks.append(mask) # insert the walk
        for i,(ngb,amask,mres,bmask) in enumerate([(my_ngb[j], and_masks[j][0], and_masks[j][1], and_masks[j][2]) for j in argsort(mdists)]):
            res.append(ngb | (amask<<4) | (mres<<10))
            bit_masks.append(bmask)

        # end point 

        quad_lists.append(mask)
        all_checks.append(len(my_ngb)) # + 2 for self mask and walk
        print(sc, ', '.join(hex(x) for x in res))
        all_masks.extend(res)
    print('// Delta to neighbour')
    print('const unsigned int delta_ngb[4] = {'+', '.join(walk_ngbs)+'};')
    print('// Hash of this')
    print('const unsigned int delta_hash[4] = {'+', '.join(hash_walk)+'};')
    hw2 = [i for sub in zip(hash_walk, walk_ngbs) for i in sub]
    print('const unsigned int hash_walk[8] = {'+', '.join(hw2)+'};')
    print('const uint16_t quad_walks[64] = {\n'+', '.join(hex(q) for q in quad_lists)+'};')
    from numpy import unique
    print('Unique walks', ','.join(hex(q) for q in unique(quad_lists)))
#    print('const unsigned char octant_masks[104] = {'+',\n'.join(', '.join(hex(x) for x in o) for o in octant_masks)+'};'
    print('Total walks', sum(all_checks))
    print('static const unsigned char ngb_start[65] = {\n'+', '.join('%d'%v for v in [0] + list(cumsum(all_checks))) + '};\n')
    print('Number of ngb_mask', len(all_masks))
    print('static const uint16_t ngb_mask[220] = {\n'+', '.join(hex(m) for m in all_masks) + '};\n')
    print('static const uint64_t ngb_bits[%d] = {\n'%len(bit_masks)+', '.join(hex(m) for m in bit_masks) + '};\n')
    print('static const uint64_t self_bits[64] = {\n'+', '.join(hex(m) for m in [find_andmask(self_block, j, nbits, dim)[2] for j in range(64)]) + '};\n')

if __name__=='__main__':
#    ngb_octants()
#    ngb_64()
    ngb2d_64()




