"""Exact E8 residual-channel certificates, extracted from the archived mixed checker.
Run through ade/run.py; see PROVENANCE.json for the source hash and changes.
"""
import os
import tempfile, shutil
from pathlib import Path
import sys
import time
import json
import math
import itertools
from collections import deque, defaultdict, Counter
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

t0 = time.time()

import e8_singleton as e8
from e7_singleton import parts, symchar

# ==============================================================================
# PART 1: E8 ORTHOGONAL-PAIR LAYER (|B| = 2, 3780 PAIRS, H_X = S6, |A_X| = 256)
# ==============================================================================
arange_240 = np.arange(240, dtype=np.uint8)
power_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 18, 20, 24, 30]

def get_power_fix_signature(P_mat):
    M = P_mat.shape[0]
    sigs = np.zeros((M, len(power_list)), dtype=np.uint16)
    powers = {1: P_mat}
    sigs[:, 0] = np.sum(P_mat == arange_240, axis=1)
    for idx_d, d in enumerate(power_list[1:], 1):
        for a in sorted(powers.keys(), reverse=True):
            if d - a in powers:
                cur = np.take_along_axis(powers[a], powers[d - a], axis=1)
                powers[d] = cur
                break
        sigs[:, idx_d] = np.sum(cur == arange_240, axis=1)
    return sigs

rep_perms = np.array([e8.word(w) for w in e8.w8], dtype=np.uint8)
rep_sigs = get_power_fix_signature(rep_perms)
sig_to_class = {tuple(int(x) for x in rep_sigs[c_idx]): c_idx for c_idx in range(112)}
assert len(sig_to_class) == 112

roots = e8.roots
lookup = e8.lookup

def canon_vec(v):
    for x in v:
        if x != 0:
            return tuple(v if x > 0 else -v)
    raise ValueError("Zero vector")

pos_mask = np.array([next(x for x in r if x != 0) > 0 for r in roots])
pos_indices = np.flatnonzero(pos_mask)
assert len(pos_indices) == 120
pos_roots = roots[pos_indices]
pos_lookup = {tuple(r): i for i, r in enumerate(pos_roots)}

refl_240 = []
refl_120 = []
for a in range(120):
    v = pos_roots[a]
    im_240 = roots - np.outer(roots @ v, v) // 4
    refl_240.append([lookup[tuple(r)] for r in im_240])
    im_120 = pos_roots - np.outer(pos_roots @ v, v) // 4
    refl_120.append([pos_lookup[canon_vec(r)] for r in im_120])

refl_240 = np.array(refl_240, dtype=np.uint8)
refl_120 = np.array(refl_120, dtype=np.int64)
gram_120 = pos_roots @ pos_roots.T

b1 = pos_lookup[(0, 0, 0, 0, 0, 0, 2, -2)]
b2 = pos_lookup[(0, 0, 0, 0, 0, 0, 2,  2)]
X0 = (min(b1, b2), max(b1, b2))
b1, b2 = X0

pairs = [(i, j) for i in range(120) for j in range(i + 1, 120) if gram_120[i, j] == 0]
assert len(pairs) == 3780
plook = {p: idx for idx, p in enumerate(pairs)}
idx_X0 = plook[X0]
simple_indices = [pos_lookup[canon_vec(r)] for r in e8.simple]

pair_mat = np.full((120, 120), -1, dtype=np.int64)
for idx, (i, j) in enumerate(pairs):
    pair_mat[i, j] = idx
    pair_mat[j, i] = idx

def ray_ref(a, v):
    v_arr = np.array(v, dtype=np.int64)
    ra = pos_roots[a]
    return canon_vec(v_arr - (int(v_arr @ ra) * ra) // 4)

def frame_ref(a, F):
    return tuple(sorted(ray_ref(a, v) for v in F))

F0 = tuple(sorted(tuple(4 * int(j == i) for j in range(8)) for i in range(6)))
frames = {idx_X0: F0}
todo = [idx_X0]
while todo:
    k = todo.pop()
    i, j = pairs[k]
    Fk = frames[k]
    for a in simple_indices:
        dest = int(pair_mat[refl_120[a, i], refl_120[a, j]])
        F_new = frame_ref(a, Fk)
        if dest in frames:
            assert frames[dest] == F_new
        else:
            frames[dest] = F_new
            todo.append(dest)
assert len(frames) == 3780

all_rays = sorted({v for F in frames.values() for v in F})
ray_to_id = {v: i for i, v in enumerate(all_rays)}
rays_arr = np.array(all_rays, dtype=np.int64)

refl_ray = np.zeros((120, len(all_rays)), dtype=np.int64)
for a in range(120):
    v = pos_roots[a]
    im_r = rays_arr - np.outer(rays_arr @ v, v) // 4
    refl_ray[a] = [ray_to_id[canon_vec(r)] for r in im_r]

frames_id = np.array([[ray_to_id[v] for v in frames[k]] for k in range(3780)], dtype=np.int64)

s6_list = list(itertools.permutations(range(6)))
s6_to_id = {p: i for i, p in enumerate(s6_list)}
s6_inv_id = np.array([s6_to_id[tuple(np.argsort(p))] for p in s6_list], dtype=np.int64)
s6_mul_id = np.array([[s6_to_id[tuple(p1[p2[k]] for k in range(6))] for p2 in s6_list] for p1 in s6_list], dtype=np.int64)
powers_6 = np.array([6**m for m in range(6)], dtype=np.int64)
code_to_s6 = np.full(6**6, -1, dtype=np.int64)
for p_id, p in enumerate(s6_list):
    code_to_s6[int(np.dot(p, powers_6))] = p_id

pairs_arr = np.array(pairs, dtype=np.int64)
pairs_i = pairs_arr[:, 0]
pairs_j = pairs_arr[:, 1]
r_dest = pair_mat[refl_120[:, pairs_i], refl_120[:, pairs_j]]

r_perm = np.zeros((120, 3780), dtype=np.int64)
for a in range(120):
    im_rays = refl_ray[a, frames_id]
    ranks = np.argsort(np.argsort(im_rays, axis=1), axis=1)
    r_perm[a] = code_to_s6[ranks @ powers_6]

out    = np.full((120, 3780), -1, dtype=np.int64)
ex     = np.zeros((120, 3780), dtype=np.int64)
e_perm = np.zeros((120, 3780), dtype=np.int64)

for a in range(120):
    for k, (i, j) in enumerate(pairs):
        if a == i or a == j:
            out[a, k] = k
            ex[a, k] = 1
            e_perm[a, k] = 0
        elif gram_120[a, i] == 0 and gram_120[a, j] == 0:
            continue
        else:
            b = i if gram_120[a, i] != 0 else j
            rd = int(r_dest[a, k])
            ed = int(r_dest[b, rd])
            out[a, k] = ed
            ex[a, k] = 0
            e_perm[a, k] = int(s6_mul_id[r_perm[b, rd], r_perm[a, k]])

B1_cols = np.zeros((3780, 76), dtype=np.int64)
B1_perms = np.zeros((3780, 76), dtype=np.int64)
B0_cols = np.zeros((3780, 1408), dtype=np.int64)
B0_perms = np.zeros((3780, 1408), dtype=np.int64)

for y_idx, (y1, y2) in enumerate(pairs):
    mid = out[y2]
    valid1 = np.flatnonzero(mid >= 0)
    mid_v = mid[valid1]
    dest = out[y1, mid_v]
    valid2 = valid1[dest >= 0]
    mid_v2 = mid[valid2]
    deg = ex[y2, valid2] + ex[y1, mid_v2]
    comb_perm = s6_mul_id[e_perm[y1, mid_v2], e_perm[y2, valid2]]
    m1 = (deg == 1)
    m0 = (deg == 0)
    B1_cols[y_idx] = valid2[m1]
    B1_perms[y_idx] = comb_perm[m1]
    B0_cols[y_idx] = valid2[m0]
    B0_perms[y_idx] = comb_perm[m0]

Z1_cols = r_dest.T.copy()
Z1_perms = np.array([[r_perm[a, Z1_cols[y_idx, a]] for a in range(120)] for y_idx in range(3780)], dtype=np.int64)

# Decompose all 3780 pairs into 12 N_W(X0)-orbits:
n_refl_roots = [a for a in range(120) if a in (b1, b2) or (gram_120[a, b1] == 0 and gram_120[a, b2] == 0)]
bridge_a = next(a for a in range(120) if gram_120[a, b1] != 0 and gram_120[a, b2] != 0)
n_gens = [(r_dest[a], r_perm[a]) for a in n_refl_roots]
d_br = np.arange(3780, dtype=np.int64)
p_br = np.zeros(3780, dtype=np.int64)
for step_a in [bridge_a, b2, b1, bridge_a]:
    p_br = s6_mul_id[r_perm[step_a, d_br], p_br]
    d_br = r_dest[step_a, d_br]
n_gens.append((d_br, p_br))

orb_id = np.full(3780, -1, dtype=np.int64)
left_perm = np.zeros(3780, dtype=np.int64)
right_perm = np.zeros(3780, dtype=np.int64)
orbit_reps = []
orbit_sizes = []

for start_k in range(3780):
    if orb_id[start_k] >= 0:
        continue
    r_idx = len(orbit_reps)
    orbit_reps.append(start_k)
    orb_id[start_k] = r_idx
    left_perm[start_k] = 0
    right_perm[start_k] = 0
    q = deque([start_k])
    cnt = 1
    while q:
        cur = q.popleft()
        lp = int(left_perm[cur])
        rp = int(right_perm[cur])
        for g_d, g_p in n_gens:
            nxt = int(g_d[cur])
            if orb_id[nxt] < 0:
                orb_id[nxt] = r_idx
                left_perm[nxt] = int(s6_mul_id[g_p[cur], lp])
                g_x0_inv = int(s6_inv_id[g_p[idx_X0]])
                right_perm[nxt] = int(s6_mul_id[rp, g_x0_inv])
                cnt += 1
                q.append(nxt)
    orbit_sizes.append(cnt)

assert len(orbit_reps) == 12 and sum(orbit_sizes) == 3780

# Method 1A: Exact N_W(X0)-orbit double-coset / character-projection trace evaluation
Z_by_deg = {1: B1_cols[idx_X0].tolist(), 0: B0_cols[idx_X0].tolist()}
gk_by_k = {}
for k in Z_by_deg[1] + Z_by_deg[0]:
    i, j = pairs[k]
    mid = int(out[b2, k])
    if ex[b2, k] == 1:
        g2 = arange_240.copy()
    else:
        c2 = i if gram_120[b2, i] != 0 else j
        g2 = refl_240[c2][refl_240[b2]]
    if ex[b1, mid] == 1:
        g1 = arange_240.copy()
    else:
        mi, mj = pairs[mid]
        c1 = mi if gram_120[b1, mi] != 0 else mj
        g1 = refl_240[c1][refl_240[b1]]
    gk_by_k[k] = g1[g2]

parts_6 = list(parts(6))
ct_to_idx = {p: idx for idx, p in enumerate(parts_6)}

def s6_cycle_type(perm6):
    seen = [False] * 6
    cyc = []
    for i in range(6):
        if not seen[i]:
            L = 0; cur = i
            while not seen[cur]:
                seen[cur] = True; cur = perm6[cur]; L += 1
            cyc.append(L)
    return tuple(sorted(cyc, reverse=True))

s6_perms_240 = []
s6_ct_idx = []
for p6 in s6_list:
    s6_ct_idx.append(ct_to_idx[s6_cycle_type(p6)])
    p8 = list(p6) + [6, 7]
    s6_perms_240.append([lookup[tuple(r)] for r in roots[:, p8]])
s6_perms_240 = np.array(s6_perms_240, dtype=np.uint8)
s6_ct_idx = np.array(s6_ct_idx, dtype=np.int64)

a_perms_240 = []
for signs in itertools.product([-1, 1], repeat=8):
    if sum(x < 0 for x in signs) % 2 != 0:
        continue
    sg = np.array(signs, dtype=np.int64)
    for swap67 in (False, True):
        out_r = roots * sg
        if swap67:
            out_r = out_r[:, [0, 1, 2, 3, 4, 5, 7, 6]]
        a_perms_240.append([lookup[tuple(r)] for r in out_r])
a_perms_240 = np.array(a_perms_240, dtype=np.uint8)

N_all = np.take_along_axis(a_perms_240[:, None, :], s6_perms_240[None, :, :], axis=2).reshape(184320, 240)
N_ct_idx = np.tile(s6_ct_idx, 256)

def compute_hist(P):
    sigs = get_power_fix_signature(P)
    uniq_sigs, inv_idx = np.unique(sigs, axis=0, return_inverse=True)
    uniq_classes = np.array([sig_to_class[tuple(int(x) for x in row)] for row in uniq_sigs], dtype=np.int64)
    e8_cids = uniq_classes[inv_idx]
    hist = np.zeros((112, 11), dtype=np.int64)
    np.add.at(hist, (e8_cids, N_ct_idx), 1)
    return hist

z1_orbits = defaultdict(list)
for k in Z_by_deg[1]:
    z1_orbits[int(orb_id[k])].append(k)
z0_orbits = defaultdict(list)
for k in Z_by_deg[0]:
    z0_orbits[int(orb_id[k])].append(k)

hist_deg1 = np.zeros((112, 11), dtype=np.int64)
for r_o, k_list in sorted(z1_orbits.items()):
    rep = k_list[0]
    gk = gk_by_k[rep]
    gk_inv = np.argsort(gk).astype(np.uint8)
    h1 = compute_hist(gk[N_all])
    h2 = compute_hist(gk_inv[N_all])
    assert np.array_equal(h1, h2)
    hist_deg1 += len(k_list) * h1

hist_deg0 = np.zeros((112, 11), dtype=np.int64)
for r_o, k_list in sorted(z0_orbits.items()):
    rep = k_list[0]
    gk = gk_by_k[rep]
    gk_inv = np.argsort(gk).astype(np.uint8)
    h1 = compute_hist(gk[N_all])
    h2 = compute_hist(gk_inv[N_all])
    assert np.array_equal(h1, h2)
    hist_deg0 += len(k_list) * h1

s6_char_table = np.array([[symchar(th, ct) for ct in parts_6] for th in parts_6], dtype=np.int64)
raw_B1_fast = e8.ch8 @ hist_deg1 @ s6_char_table.T
raw_B0_fast = e8.ch8 @ hist_deg0 @ s6_char_table.T
raw_B1_exact = e8.ch8.astype(object) @ hist_deg1.astype(object) @ s6_char_table.T.astype(object)
raw_B0_exact = e8.ch8.astype(object) @ hist_deg0.astype(object) @ s6_char_table.T.astype(object)
assert np.array_equal(raw_B1_fast.astype(object), raw_B1_exact)
assert np.array_equal(raw_B0_fast.astype(object), raw_B0_exact)
print('Exact object-integer trace products match int64:', max(abs(int(x)) for x in raw_B1_exact.flat), max(abs(int(x)) for x in raw_B0_exact.flat), flush=True)
raw_B1 = raw_B1_fast
raw_B0 = raw_B0_fast
assert np.all(raw_B1 % 184320 == 0) and np.all(raw_B0 % 184320 == 0)
tr_B1_iso = raw_B1 // 184320
tr_B0_iso = raw_B0 // 184320

with open('e8_pair_branching_results.json') as f:
    B_pair = np.array(json.load(f)['branching'], dtype=np.int64)

assert np.all(tr_B1_iso[B_pair == 0] == 0)
assert np.all(tr_B0_iso[B_pair == 0] == 0)
assert np.all((e8.ch8[:, 0] @ tr_B1_iso) == 0)
assert np.all((e8.ch8[:, 0] @ tr_B0_iso) == 0)

# Method 1B: Direct 12-orbit compressed Z_1-Lagrange projection & integral Specht eigenvectors
def get_syts(part):
    cells = [(r, c) for r, row_len in enumerate(part) for c in range(row_len)]
    syts = []
    def backtrack(placed):
        if len(placed) == 6:
            syts.append(tuple(placed))
            return
        for r, c in cells:
            if (r, c) not in placed and (r == 0 or (r - 1, c) in placed) and (c == 0 or (r, c - 1) in placed):
                backtrack(placed + [(r, c)])
    backtrack([])
    return syts

s6_words = {0: ()}
q_s6 = deque([0])
adj_gens = []
for k in range(5):
    g = list(range(6)); g[k], g[k+1] = g[k+1], g[k]
    adj_gens.append(s6_to_id[tuple(g)])
while q_s6:
    u = q_s6.popleft()
    w = s6_words[u]
    for k, g_id in enumerate(adj_gens):
        v = int(s6_mul_id[u, g_id])
        if v not in s6_words:
            s6_words[v] = w + (k,)
            q_s6.append(v)

def build_integral_specht_720(part):
    syts = get_syts(part)
    d = len(syts)
    syt_idx = {t: i for i, t in enumerate(syts)}
    gen_mats = []
    for k in range(5):
        M = sp.zeros(d, d)
        for i, t in enumerate(syts):
            r1, c1 = t[k]; r2, c2 = t[k+1]
            axial = (c2 - r2) - (c1 - r1)
            t_swap = list(t); t_swap[k], t_swap[k+1] = t_swap[k+1], t_swap[k]; t_swap = tuple(t_swap)
            if t_swap not in syt_idx:
                M[i, i] = sp.Rational(1, axial)
            else:
                j = syt_idx[t_swap]
                if i < j:
                    M[i, i] = sp.Rational(1, axial)
                    M[j, i] = sp.Integer(1)
                    M[i, j] = sp.Integer(1) - sp.Rational(1, axial * axial)
                    M[j, j] = -sp.Rational(1, axial)
        gen_mats.append(M)
    mats = [None] * 720
    mats[0] = sp.eye(d)
    for p_id in range(1, 720):
        w = s6_words[p_id]
        prev_id = int(s6_mul_id[p_id, adj_gens[w[-1]]])
        mats[p_id] = mats[prev_id] * gen_mats[w[-1]]
    cols = [mats[p][:, 0] for p in range(720)]
    L = 1
    for c in cols:
        for x in c:
            L = math.lcm(L, int(x.q))
    int_cols = [[int(c[r] * L) for c in cols] for r in range(d)]
    H = hermite_normal_form(sp.Matrix(int_cols))
    B_rat = H / sp.Integer(L)
    B_inv = B_rat.inv()
    int_gen_mats = []
    for k in range(5):
        G_int = B_inv * gen_mats[k] * B_rat
        assert all(x.is_integer for x in G_int)
        int_gen_mats.append(np.array([[int(G_int[r, c]) for c in range(d)] for r in range(d)], dtype=np.int64))
    int_mats = np.zeros((720, d, d), dtype=np.int64)
    int_mats[0] = np.eye(d, dtype=np.int64)
    for p_id in range(1, 720):
        w = s6_words[p_id]
        prev_id = int(s6_mul_id[p_id, adj_gens[w[-1]]])
        exact_product = int_mats[prev_id].astype(object) @ int_gen_mats[w[-1]].astype(object)
        assert all(-(2**63) <= int(x) < 2**63 for x in exact_product.flat)
        int_mats[p_id] = exact_product
    return int_mats, int_gen_mats

def compile_compressed_op(cols_t, perms_t):
    rep_rows = np.array(orbit_reps, dtype=np.int64)
    c_sub = cols_t[rep_rows]
    p_sub = perms_t[rep_rows]
    o_sub = orb_id[c_sub]
    l_sub = s6_mul_id[p_sub, left_perm[c_sub]]
    r_sub = right_perm[c_sub]
    return o_sub, l_sub, r_sub

Z1_comp = compile_compressed_op(Z1_cols, Z1_perms)
B1_comp = compile_compressed_op(B1_cols, B1_perms)
B0_comp = compile_compressed_op(B0_cols, B0_perms)

def apply_compressed_obj(comp_op, S_mats_obj, V_rep_obj):
    o_sub, l_sub, r_sub = comp_op
    d = S_mats_obj.shape[1]
    out_rep = np.zeros((12, d, d), dtype=object)
    for r in range(12):
        acc = np.zeros((d, d), dtype=object)
        for r_id in np.unique(r_sub[r]):
            mask = (r_sub[r] == r_id)
            o_m = o_sub[r, mask]
            l_m = l_sub[r, mask]
            lhs = np.sum(S_mats_obj[l_m] @ V_rep_obj[o_m], axis=0)
            acc += lhs @ S_mats_obj[r_id]
        out_rep[r] = acc
    return out_rep

def gcd_array_obj(arr):
    g = 0
    for x in arr.flat:
        g = math.gcd(g, abs(int(x)))
        if g == 1:
            return 1
    return g

r_X0 = int(orb_id[idx_X0])
eigenvector_certificates = []

for part, targets in [
    ((4, 2), [(18, -3, [(2240, 10)])]),
    ((2, 2, 1, 1), [(-18, 5, [(2240, 28)])]),
    ((3, 2, 1), [(8, -3, [(2835, 14), (6075, 14)]), (-8, 5, [(2835, 22), (6075, 22)])]),
]:
    j_part = parts_6.index(part)
    spec_Z1 = sorted({int(e8.c8[i]) for i in range(112) if B_pair[i, j_part] > 0})
    S_mats, int_gen_mats = build_integral_specht_720(part)
    S_mats_obj = S_mats.astype(object)
    d = S_mats.shape[1]
    I_rep = np.zeros((12, d, d), dtype=object)
    I_rep[r_X0] = np.eye(d, dtype=object)

    P = I_rep.copy()
    for c in spec_Z1:
        P = apply_compressed_obj(Z1_comp, S_mats_obj, P) - c * P
        g = gcd_array_obj(P)
        if g > 1:
            P //= g
    assert all(int(x) == 0 for x in P.flat)
    print('Annihilating polynomial verified for', part, flush=True)

    for c0, delta_target, w_labels in targets:
        V_rep = I_rep.copy()
        for c in spec_Z1:
            if c == c0:
                continue
            ZV_rep = apply_compressed_obj(Z1_comp, S_mats_obj, V_rep)
            V_rep = ZV_rep - c * V_rep
            g = gcd_array_obj(V_rep)
            if g > 1:
                V_rep //= g
        kappa = int(V_rep[r_X0, 0, 0])
        assert kappa != 0
        assert np.array_equal(V_rep[r_X0], kappa * np.eye(d, dtype=object))

        Z1_V = apply_compressed_obj(Z1_comp, S_mats_obj, V_rep)
        assert np.array_equal(Z1_V, c0 * V_rep)
        B1_V = apply_compressed_obj(B1_comp, S_mats_obj, V_rep)
        B0_V = apply_compressed_obj(B0_comp, S_mats_obj, V_rep)
        assert B1_V[r_X0, 0, 0] % kappa == 0 and B0_V[r_X0, 0, 0] % kappa == 0
        a_val = int(B1_V[r_X0, 0, 0] // kappa)
        b_val = int(B0_V[r_X0, 0, 0] // kappa)
        assert np.array_equal(B1_V, a_val * V_rep)
        assert np.array_equal(B0_V, b_val * V_rep)

        # Cross-check against Method 1A character-projection traces:
        for wl in w_labels:
            i_w = next(idx for idx, p in enumerate(e8.pa8) if tuple(p) == wl)
            assert int(tr_B1_iso[i_w, j_part]) == a_val
            assert int(tr_B0_iso[i_w, j_part]) == b_val

        D_sym = sp.Symbol('delta')
        q_fact = str(sp.factor(D_sym**2 + a_val * D_sym + b_val))
        q_eval = delta_target**2 + a_val * delta_target + b_val
        assert q_eval != 0
        eigenvector_certificates.append({
            "delta_target": delta_target,
            "S6_partition": list(part),
            "S6_dim": d,
            "W_E8_labels": [list(wl) for wl in w_labels],
            "Z1_eigenvalue_c0": c0,
            "Spec_Z1_on_module": spec_Z1,
            "base_fiber_scalar_kappa": kappa,
            "max_abs_entry_V_rep": max(abs(int(x)) for x in V_rep.flat),
            "a_coeff": a_val,
            "b_coeff": b_val,
            "quadratic_factored": q_fact,
            "value_at_delta_target": q_eval,
            "integral_specht_generators": [g.tolist() for g in int_gen_mats],
            "compressed_eigenvector_12_orbits_col0": [[int(x) for x in V_rep[r, :, 0]] for r in range(12)],
        })

# Collect all 115 multiplicity-1 quadratics across all 11 S6 partitions:
all_mult1_table = []
all_mult1_roots = set()
D_sym = sp.Symbol('delta')
for j, th in enumerate(parts_6):
    for i in range(112):
        if B_pair[i, j] == 1:
            a_val = int(tr_B1_iso[i, j])
            b_val = int(tr_B0_iso[i, j])
            disc = a_val * a_val - 4 * b_val
            sr = int(round(math.isqrt(disc))) if disc >= 0 else -1
            assert disc >= 0 and sr * sr == disc and (-a_val + sr) % 2 == 0
            r1, r2 = sorted(((-a_val - sr) // 2, (-a_val + sr) // 2))
            for r in (r1, r2):
                if r != 0:
                    all_mult1_roots.add(r)
            all_mult1_table.append({
                "S6_partition": list(th),
                "W_E8_label": list(e8.pa8[i]),
                "c_W": int(e8.c8[i]),
                "a": a_val,
                "b": b_val,
                "roots": [r1, r2],
            })

S8_set = set(e8.result['singleton_bad_set'])
assert len(all_mult1_table) == 115
assert all_mult1_roots <= S8_set

# ==============================================================================
payload = {
    "status": "all exact integer and rational assertions passed (zero floating-point tolerances)",
    "elapsed_seconds": round(time.time() - t0, 2),
    "E8_residual_certificates": {
        "nonzero_exceptional_set_size": 40,
        "nonzero_exceptional_set_S8": sorted(S8_set),
        "decided_candidates_minus3_and_5": "all 6 residual channels nonzero; the paper supplies the reduction",
        "N_W_X0_order": 184320,
        "A_X0_order": 256,
        "H_X0": "S_6 (order 720)",
        "pair_orbit_size": 3780,
        "N_W_X0_orbit_sizes_on_3780_pairs": orbit_sizes,
        "Z1_degree1_block_count": 76,
        "Z1_orbit_sizes": [len(v) for _, v in sorted(z1_orbits.items())],
        "Z0_degree0_block_count": 1408,
        "Z0_orbit_sizes": [len(v) for _, v in sorted(z0_orbits.items())],
        "six_residual_channel_certificates": eigenvector_certificates,
        "total_multiplicity_one_constituents": len(all_mult1_table),
        "all_115_multiplicity_one_quadratics_split_over_Z": True,
        "union_of_all_115_multiplicity_one_nonzero_roots": sorted(all_mult1_roots),
        "multiplicity_one_roots_outside_S8": sorted(all_mult1_roots - S8_set),
        "multiplicity_one_table_115": all_mult1_table,
        "joint_histogram_deg1_112x11": hist_deg1.tolist(),
        "joint_histogram_deg0_112x11": hist_deg0.tolist(),
        "isotypic_trace_B1_112x11": tr_B1_iso.tolist(),
        "isotypic_trace_B0_112x11": tr_B0_iso.tolist(),
    },
}
with open("annihilator_verification.json","w") as f:
    json.dump(payload,f,indent=2)
print("PASS exact E8 residual certificates",payload["elapsed_seconds"])
