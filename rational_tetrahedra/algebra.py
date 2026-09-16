"""Exact Laurent and cyclotomic arithmetic used by the audit.

Only the Python standard library is used.  Laurent polynomials are sparse
dictionaries whose keys are exponent tuples.  Ordinary polynomials are lists
of integer coefficients in increasing degree order.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import gcd
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

Exponent = Tuple[int, int, int, int, int, int]
Laurent = Dict[Exponent, Fraction]

# The order used in the paper pairs opposite edges.
EDGES: Tuple[Tuple[int, int], ...] = (
    (0, 1),  # 12
    (2, 3),  # 34
    (0, 2),  # 13
    (1, 3),  # 24
    (0, 3),  # 14
    (1, 2),  # 23
)
EDGE_INDEX = {tuple(sorted(edge)): i for i, edge in enumerate(EDGES)}
ZERO_EXP: Exponent = (0, 0, 0, 0, 0, 0)


def _add_term(poly: Laurent, exponent: Exponent, coefficient: Fraction) -> None:
    value = poly.get(exponent, Fraction(0)) + coefficient
    if value:
        poly[exponent] = value
    elif exponent in poly:
        del poly[exponent]


def _mul(a: Mapping[Exponent, Fraction], b: Mapping[Exponent, Fraction]) -> Laurent:
    out: Laurent = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            exponent = tuple(x + y for x, y in zip(ea, eb))
            _add_term(out, exponent, ca * cb)
    return out


def _entry(i: int, j: int) -> Laurent:
    if i == j:
        return {ZERO_EXP: Fraction(1)}
    k = EDGE_INDEX[tuple(sorted((i, j)))]
    positive = [0] * 6
    negative = [0] * 6
    positive[k] = 1
    negative[k] = -1
    return {
        tuple(positive): Fraction(1, 2),
        tuple(negative): Fraction(1, 2),
    }


def _permutation_sign(p: Sequence[int]) -> int:
    inversions = sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))
    return -1 if inversions % 2 else 1


@lru_cache(maxsize=1)
def gram_laurent_polynomial() -> Laurent:
    """Return det(cos(theta_ij)) after z_ij = exp(i theta_ij).

    The result has 105 Laurent monomials.  Multiplication by 16 makes all
    coefficients integral and gives constant coefficient -20, matching
    equation (2) in the paper.
    """

    result: Laurent = {}
    for p in permutations(range(4)):
        term: Laurent = {ZERO_EXP: Fraction(_permutation_sign(p))}
        for i, j in enumerate(p):
            term = _mul(term, _entry(i, j))
        for exponent, coefficient in term.items():
            _add_term(result, exponent, coefficient)
    return result


def _lcm(a: int, b: int) -> int:
    return abs(a // gcd(a, b) * b) if a and b else 0


@lru_cache(maxsize=1)
def integral_gram_terms() -> Tuple[int, Tuple[Tuple[Exponent, int], ...]]:
    """Return the smallest positive scale and integral Laurent terms."""

    poly = gram_laurent_polynomial()
    scale = 1
    for coefficient in poly.values():
        scale = _lcm(scale, coefficient.denominator)
    terms = tuple(sorted((e, int(c * scale)) for e, c in poly.items()))
    return scale, terms


def _trim(poly: List[int]) -> List[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def _div_exact(numerator: Sequence[int], denominator: Sequence[int]) -> List[int]:
    """Exact division in Z[x] for a monic denominator."""

    num = list(numerator)
    den = _trim(list(denominator))
    if den[-1] != 1:
        raise ValueError("denominator must be monic")
    if len(num) < len(den):
        raise ValueError("non-exact polynomial division")
    quotient = [0] * (len(num) - len(den) + 1)
    while len(num) >= len(den):
        coefficient = num[-1]
        shift = len(num) - len(den)
        quotient[shift] = coefficient
        for i, value in enumerate(den):
            num[i + shift] -= coefficient * value
        _trim(num)
    if any(num):
        raise ValueError("non-exact polynomial division")
    return _trim(quotient)


def _proper_divisors(n: int) -> Iterable[int]:
    return (d for d in range(1, n) if n % d == 0)


@lru_cache(maxsize=None)
def cyclotomic_polynomial(n: int) -> Tuple[int, ...]:
    """Return Phi_n with coefficients in increasing degree order."""

    if n < 1:
        raise ValueError("cyclotomic index must be positive")
    poly = [-1] + [0] * (n - 1) + [1]
    for d in _proper_divisors(n):
        poly = _div_exact(poly, cyclotomic_polynomial(d))
    return tuple(poly)


def polynomial_remainder(poly: Sequence[int], monic_modulus: Sequence[int]) -> Tuple[int, ...]:
    """Return the remainder of an integer polynomial modulo a monic one."""

    value = _trim(list(poly) or [0])
    modulus = _trim(list(monic_modulus))
    if modulus[-1] != 1:
        raise ValueError("modulus must be monic")
    while len(value) >= len(modulus):
        coefficient = value[-1]
        shift = len(value) - len(modulus)
        for i, m in enumerate(modulus):
            value[i + shift] -= coefficient * m
        _trim(value)
    return tuple(value)


def exact_torsion_zero(
    level: int,
    numerators: Sequence[int],
    angle_kind: str = "dihedral",
) -> bool:
    """Certify the Gram equation for angles pi*numerator/level.

    The Laurent polynomial is written for angles between face-normal lines.
    For an interior dihedral angle alpha, outward normals meet at pi-alpha,
    so their phase is q**(level-m), not q**m.  ``angle_kind`` makes this
    convention explicit and may be either ``"dihedral"`` or ``"line"``.
    """

    if level < 2 or len(numerators) != 6:
        raise ValueError("expected a level >= 2 and six angle numerators")
    if angle_kind not in {"dihedral", "line"}:
        raise ValueError("angle_kind must be 'dihedral' or 'line'")
    order = 2 * level
    phase_powers = (
        tuple(level - m for m in numerators)
        if angle_kind == "dihedral"
        else tuple(numerators)
    )
    _, terms = integral_gram_terms()
    accumulated = [0] * order
    for exponent, coefficient in terms:
        power = sum(e * m for e, m in zip(exponent, phase_powers)) % order
        accumulated[power] += coefficient
    remainder = polynomial_remainder(accumulated, cyclotomic_polynomial(order))
    return not any(remainder)


def _parametric_identity(
    root_order: int,
    constant_powers: Sequence[int],
    parameter_powers: Sequence[int],
) -> bool:
    """Verify a Laurent identity with z_j=q**a_j*u**b_j exactly."""

    if len(constant_powers) != 6 or len(parameter_powers) != 6:
        raise ValueError("six substitutions are required")
    _, terms = integral_gram_terms()
    by_parameter: Dict[int, List[int]] = {}
    for exponent, coefficient in terms:
        upower = sum(e * b for e, b in zip(exponent, parameter_powers))
        qpower = sum(e * a for e, a in zip(exponent, constant_powers)) % root_order
        coeffs = by_parameter.setdefault(upower, [0] * root_order)
        coeffs[qpower] += coefficient
    modulus = cyclotomic_polynomial(root_order)
    return all(
        not any(polynomial_remainder(coeffs, modulus))
        for coeffs in by_parameter.values()
    )


def verify_parametric_families() -> Mapping[str, bool]:
    """Exactly certify the two families in Theorem 1.8.

    Put q=exp(i*pi/6), a primitive 12th root, and u=exp(i*x).
    Each angle phase is represented as q**a * u**b.
    """

    # Convert each interior dihedral phase q**a*u**b to the outward-normal
    # phase -q**(-a)*u**(-b) = q**(6-a)*u**(-b).
    family_1 = _parametric_identity(
        12,
        constant_powers=(3, 3, 0, 4, 6, 6),
        parameter_powers=(0, 0, 2, 0, -1, -1),
    )
    family_2 = _parametric_identity(
        12,
        constant_powers=(1, 5, 2, 2, 6, 6),
        parameter_powers=(1, -1, 1, 1, -1, -1),
    )
    return {"family_1": family_1, "family_2": family_2}
