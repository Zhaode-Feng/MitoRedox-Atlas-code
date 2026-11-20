"""Original BLOSUM62 and Gotoh optimal-centre mapping algorithm."""

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62_TEXT = """
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4
"""
BLOSUM = {}
for _line in BLOSUM62_TEXT.strip().splitlines():
    _values = _line.split()
    for _aa, _score in zip(AA_ORDER, _values[1:]):
        BLOSUM[_values[0], _aa] = int(_score)
STANDARD_AA = frozenset(AA_ORDER)

class ContractError(RuntimeError):
    pass


class ResourceCap(ContractError):
    pass

def window13(sequence, position):
    if not 1 <= position <= len(sequence):
        raise ContractError("WINDOW_POSITION_OUT_OF_RANGE")
    start, stop = max(1, position - 6), min(len(sequence), position + 6)
    return {"wt_human_window": sequence[start - 1:stop], "window_start": start,
            "window_stop": stop, "window_length": stop - start + 1,
            "window_complete_13aa": stop - start + 1 == 13}


def alignment_primary(identical, nongap, n, m, assignments):
    return (len(assignments) == 1 and 0 not in assignments and nongap > 0
            and 100 * identical >= 50 * nongap
            and 100 * nongap >= 80 * n and 100 * nongap >= 80 * m)


def global_alignment(mouse, human, centers, check=lambda: None, pair_cap=2_000_000):
    """Gotoh -11/-1. Store every optimal predecessor; reverse all optimal paths.
    Reverse reachability is exactly the union-of-center-assignments DP contract.
    X consumes mouse/gap in human; Y consumes human/gap in mouse.
    X opens only from M, Y opens only from M; opposite adjacent gaps are disallowed.
    Deterministic traceback chooses M > X > Y, without changing ambiguity.
    """
    if not set(mouse + human).issubset(STANDARD_AA):
        raise ContractError("NONSTANDARD_SEQUENCE_ALIGNMENT_DEFERRED")
    n, m = len(mouse), len(human)
    if n * m > pair_cap:
        raise ResourceCap("DEFERRED_RESOURCE_CAP")
    if any(not 1 <= c <= n for c in centers):
        raise ContractError("CENTER_OUT_OF_RANGE")
    width, cells, neg = m + 1, (n + 1) * (m + 1), -10**9
    masks = [bytearray(cells) for _ in range(3)]
    pm, px, py = [neg] * width, [neg] * width, [neg] * width
    pm[0] = 0
    for j in range(1, width):
        py[j] = -11 - (j - 1)
        masks[2][j] = 4 if j > 1 else 1
    for i in range(1, n + 1):
        cm, cx, cy = [neg] * width, [neg] * width, [neg] * width
        cx[0] = -11 - (i - 1)
        masks[1][i * width] = 2 if i > 1 else 1
        base = i * width
        a = mouse[i - 1]
        for j in range(1, width):
            best = max(pm[j - 1], px[j - 1], py[j - 1])
            cm[j] = best + BLOSUM[a, human[j - 1]]
            masks[0][base + j] = ((pm[j - 1] == best) + 2 * (px[j - 1] == best)
                                  + 4 * (py[j - 1] == best))
            opening, extending = pm[j] - 11, px[j] - 1
            cx[j] = max(opening, extending)
            masks[1][base + j] = (opening == cx[j]) + 2 * (extending == cx[j])
            opening, extending = cm[j - 1] - 11, cy[j - 1] - 1
            cy[j] = max(opening, extending)
            masks[2][base + j] = (opening == cy[j]) + 4 * (extending == cy[j])
        pm, px, py = cm, cx, cy
        if i % 16 == 0:
            check()
    best = max(pm[m], px[m], py[m])
    terminal = (pm[m] == best) + 2 * (px[m] == best) + 4 * (py[m] == best)
    reachable = bytearray(cells)
    reachable[-1] = terminal
    assignments = {c: set() for c in centers}
    tied = bool(terminal & (terminal - 1))
    for i in range(n, -1, -1):
        for j in range(m, -1, -1):
            cell = i * width + j
            states = reachable[cell]
            if not states:
                continue
            for state, bit in enumerate((1, 2, 4)):
                if not states & bit or i == j == 0:
                    continue
                predecessor = masks[state][cell]
                tied |= bool(predecessor & (predecessor - 1))
                if state == 0:
                    if i in assignments:
                        assignments[i].add(j)
                    parent = cell - width - 1
                elif state == 1:
                    if i in assignments:
                        assignments[i].add(0)
                    parent = cell - width
                else:
                    parent = cell - 1
                if parent < 0:
                    raise ContractError("INVALID_OPTIMAL_TRACE")
                reachable[parent] |= predecessor
        if i % 16 == 0:
            check()
    state = next(s for s, bit in enumerate((1, 2, 4)) if terminal & bit)
    i, j = n, m
    left, right, mapping = [], [], {}
    identical = nongap = 0
    while i or j:
        cell = i * width + j
        pred = masks[state][cell]
        if state == 0:
            left.append(mouse[i - 1]); right.append(human[j - 1])
            mapping[i] = j
            nongap += 1
            identical += mouse[i - 1] == human[j - 1]
            i -= 1; j -= 1
        elif state == 1:
            left.append(mouse[i - 1]); right.append("-"); mapping[i] = None; i -= 1
        else:
            left.append("-"); right.append(human[j - 1]); j -= 1
        state = next(s for s, bit in enumerate((1, 2, 4)) if pred & bit)
    return {"score": best, "identical_pairs": identical, "nongap_pairs": nongap,
            "identity_percent": 100 * identical / nongap if nongap else 0,
            "mouse_coverage_percent": 100 * nongap / n if n else 0,
            "human_coverage_percent": 100 * nongap / m if m else 0,
            "center_assignments": {str(c): sorted(values) for c, values in assignments.items()},
            "global_optimal_tie": tied, "deterministic_mapping": mapping,
            "aligned_mouse": "".join(reversed(left)), "aligned_human": "".join(reversed(right)),
            "dp_cells": n * m, "allocated_matrix_cells": cells}
