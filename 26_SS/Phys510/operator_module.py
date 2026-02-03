import sympy as sp

def split_scalar_operator(term):
    coeff = sp.S.One
    ops = []

    for f in sp.Mul.make_args(term):
        if f.is_commutative:
            coeff *= f
        else:
            ops.append(f)
    op_prod = sp.Mul(*ops) if ops else sp.S.One

    return sp.simplify(coeff), op_prod

def canonical_commute(op_prod, commuting_pairs = None, anticommuting_pairs = None):
    if op_prod is sp.S.One:
        return op_prod, 1

    ops = list(sp.Mul.make_args(op_prod))
    sign = 1
    changed = True
    while changed:
        changed = False
        for i in range(len(ops) - 1):
            a, b = ops[i], ops[i+1]
            key = frozenset([a, b])
            if commuting_pairs and key in commuting_pairs:
                    if sp.default_sort_key(a) > sp.default_sort_key(b): # the ordering that is auto-picked
                        ops[i], ops[i+1] = b, a
                        changed = True
            elif anticommuting_pairs and key in anticommuting_pairs:
                    if sp.default_sort_key(a) > sp.default_sort_key(b): # the ordering that is auto-picked
                        ops[i], ops[i+1] = b, a
                        sign *= -1
                        changed = True
    return sp.Mul(*ops), sign
    

def collect_operators(expr, collect_coeffs = False, commuting_pairs = None, anticommuting_pairs = None):
    expr = sp.expand(expr)
    terms = sp.Add.make_args(expr)
    collected = {}

    for t in terms:
        sign = 1
        coeff, op = split_scalar_operator(t)
        if commuting_pairs or anticommuting_pairs:
            op, sign = canonical_commute(op, commuting_pairs, anticommuting_pairs)
        collected[op] = collected.get(op, 0) + (coeff * sign)
    if not collect_coeffs:
        return collected

    # invert dict to collect obvious coefficients
    coeff_collected = {}
    for op, coeff in collected.items():
        coeff = sp.simplify(coeff)
        coeff_collected.setdefault(coeff, 0)
        coeff_collected[coeff] += op

    return coeff_collected

def replace_operator(expr, old, new):

    result = 0
    for term in sp.Add.make_args(sp.expand(expr)):
        coeff, op = split_scalar_operator(term)
        if sp.simplify(op-old)==0:
            result += coeff*new            
        else:
            result += coeff*op
    return result
