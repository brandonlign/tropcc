"""Factor explicit mass conditions; this does not extract the theorem's full q.

Positive-zero witnesses describe certificate failures, not infinite central
configurations. Unclassified factors are retained rather than assumed absent.
"""
from pathlib import Path
from collections import Counter
import argparse, hashlib, json
import flint
import sympy as sp
import verify as V

ROOT = Path(__file__).resolve().parents[1]
BATCHES = ['cell_binomial', 'equality', 'shifted', 'reduction', 'target16',
           'target301', 'pair', 'transfer', 'pair2', 'transfer2', 'target247',
           'target279', 'target897', 'transfer897']

def positive_status(f):
    coefficients = list(f.to_dict().values())
    if all(c > 0 for c in coefficients) or all(c < 0 for c in coefficients):
        return {'status': 'nonvanishing_same_sign_coefficients'}
    grid = [flint.fmpq(1, 8), flint.fmpq(1, 4), flint.fmpq(1, 2)]
    grid += [flint.fmpq(x) for x in [1, 2, 3, 4, 7, 8, 16, 32, 64]]
    signs = {}
    for a in grid:
        for b in grid:
            value = f(a, b)
            if value == 0:
                return {'status': 'positive_zero', 'witness': [str(a), str(b)]}
            signs.setdefault(1 if value > 0 else -1,
                             [str(a), str(b), str(value)])
            if len(signs) == 2:
                # The positive quadrant is convex; the polynomial is continuous.
                return {'status': 'positive_zero_by_continuity',
                        'opposite_sign_endpoints': list(signs.values())}
    return {'status': 'undetermined'}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'logs/mass_factor_inventory.json')
    args = parser.parse_args()
    ctx = flint.fmpq_mpoly_ctx.get(('a', 'b'), 'lex')
    expressions = {}; inputs = {}; terms = 0
    for folder in ['cert', 'dcert', 'gcert', 'hcert', 'ecert', 'split']:
        for path in sorted((ROOT / 'data/exact' / folder).glob('*')):
            if path.suffix not in ('.cof', '.split'):
                continue
            relative = str(path.relative_to(ROOT))
            inputs[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
            with path.open() as stream:
                for line in stream:
                    if '|' in line:
                        parts = line.strip().split('|'); assert len(parts) == 3
                        expressions.setdefault(parts[1], set()).add(relative)
                        terms += 1
    factors = {}
    for i, (text, files) in enumerate(expressions.items()):
        f = V.parse_poly(text, ctx, ['a', 'b']); assert f != 0
        unit, decomposition = f.factor()
        reconstructed = ctx.constant(unit)
        for g, power in decomposition:
            reconstructed *= g ** power
            record = factors.setdefault(str(g), {'polynomial': str(g),
                                                  'files': set()})
            record['files'].update(files)
        assert reconstructed == f, 'factorization reconstruction mismatch'
        if i % 5000 == 0:
            print('Factored', i, 'of', len(expressions), flush=True)
    explicit = []
    for text, record in sorted(factors.items()):
        f = V.parse_poly(text, ctx, ['a', 'b'])
        explicit.append(dict(record, files=sorted(record['files']),
                             total_degree=int(max(sum(m) for m in f.to_dict())),
                             positive=positive_status(f)))
    a, b, c = sp.symbols('a b c'); companion = {}; identities = 0
    for batch in BATCHES:
        path = ROOT / 'companion/results' / (batch + '_certificates.json')
        inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        for record in json.loads(path.read_text())['certificates']:
            identities += 1
            companion.setdefault(record['parameter_factor'], set()).add(batch)
    companion_factors = {}
    for text, batches in companion.items():
        f = sp.Poly(sp.sympify(text).subs(c, 1), a, b, domain=sp.QQ)
        assert not f.is_zero
        unit, decomposition = sp.factor_list(f); reconstructed = sp.Poly(unit, a, b)
        for g, power in decomposition:
            reconstructed *= g ** power
            _, primitive = g.clear_denoms()
            _, primitive = primitive.primitive()
            polynomial = ctx.from_dict({m: int(v) for m, v in primitive.terms()})
            companion_factors.setdefault(str(polynomial), set()).update(batches)
        assert reconstructed == f
    companion_records = [dict(polynomial=text, batches=sorted(batches),
                              positive=positive_status(V.parse_poly(text, ctx, ['a', 'b'])))
                         for text, batches in sorted(companion_factors.items())]
    candidates = []
    for aa, bb in [(3, 7), (2, 3), (2, 5), (3, 5), (5, 11), (7, 13), (11, 17)]:
        candidates.append({'mass_pair': [aa, bb], 'certifies_instance': False,
                           'cofactor_zero_factors': [r['polynomial'] for r in explicit
                               if V.parse_poly(r['polynomial'], ctx, ['a', 'b'])(aa, bb) == 0],
                           'companion_zero_factors': [r['polynomial'] for r in companion_records
                               if V.parse_poly(r['polynomial'], ctx, ['a', 'b'])(aa, bb) == 0]})
    report = {'full_exceptional_polynomial_extracted': False,
              'cofactor_files': len(inputs) - len(BATCHES), 'cofactor_terms': terms,
              'distinct_denominators': len(expressions),
              'explicit_factor_count': len(explicit),
              'positive_counts': dict(Counter(r['positive']['status'] for r in explicit)),
              'explicit_factors': explicit, 'companion_identities': identities,
              'companion_distinct_expressions': len(companion),
              'companion_factors_at_c_1': companion_records,
              'rational_candidates': candidates, 'input_sha256': inputs,
              'scope': 'Explicit identity denominators and companion factors only. '
                       'Does not extract existential modular/elimination spreading '
                       'conditions. Positive zeros are failures of recorded identities, '
                       'not evidence of infinite configurations. Multiple proof routes '
                       'may remove some of these conditions.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print('PASS factor reconstruction and exact positive witnesses:', report['positive_counts'])
    print('Companion irreducible factors:', len(companion_records))

if __name__ == '__main__':
    main()
