# Infinitely many general rank sequences

This note answers the fixed-rank finiteness conjecture at arXiv:2501.14643v1, Section 1, page 2. At every initial rank `r >= 3`, positive integer sequences realize infinitely many general rank sequences. A second family gives infinitely many particular rank-three rank sequences under one fixed recurrence.

## Mathematical contents

- Theorem 1 gives the fixed-recurrence family `s_p(n) = p*8^n - 4^n - 2^n`, for `p >= 3`, and its coefficient-support formula.
- Theorem 5 gives `u_(r,K)(n) = 1 + 2^n + 2^(K*n) + sum(q_i^n)` for `r-3` distinct odd primes, and its exact power rank `binomial(M+r-1,r-1) - binomial(M-K+r-1,r-1)`, the second term being zero for `M < K`.
- Corollaries give arbitrarily long shared prefixes, the arbitrary-base version of `s_p`, and exact homogeneous relation groups.
- The coverage table credits the posted polynomiality argument and distinguishes the particular-class assertion of Conjecture 5.3 from the general-class result.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a current LaTeX distribution:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is contained in `note.tex`; no separate bibliography file is needed. Tectonic may retrieve its standard TeX bundle on its first compilation.

## Exact finite checks

Python 3.10 or later and the standard library suffice:

```sh
python3 verify.py > verification.json
```

The checker compares convolution with multinomial coefficients for every `3 <= p <= 30` and `1 <= M <= 65`. It checks forced cancellation at `M=2*p+1`, 4,200 exponential identities, positive trinomial coefficients, five exact Hankel determinants, and annihilating recurrences. For `3 <= r <= 6`, `2 <= K <= 12`, and `1 <= M <= 12`, it constructs positive-coefficient product expansions and compares their support sizes with the general-family formula. It also checks 192 additional sequence identities and 44 first-departure comparisons. All infinite-family conclusions have written proofs in `note.tex`.

The expected SHA-256 of the generated output is recorded in `checks.sha256`.

## Formal statement and written scope

The [pinned Lean file](https://github.com/the-omega-institute/trureturing/blob/b0f7c01b01511ea00428c1b6edad5dd2b372eaf8/D5/S1/Recurrence/RowlandSistosBarronRankSequences.lean) contains `D5.S1.Recurrence.RowlandSistosBarronRankSequences.result : ¬ claim`, where `claim` is fixed-rank finiteness for integer sequences with rational tail recurrences and all positive exponents. Its proof uses the fixed rank-three family `s_p` and finite fibres. Positivity, the common recurrence, the complex-coefficient lower bound, the general-family formula, and all corollaries have written proofs in the note. The finite checks have exactly the ranges stated above.

## Files and license

The source archive contains `note.tex`, `metadata.json`, `verify.py`, `README.md`, and `checks.sha256`. The PDF and generated checker output are supplied separately.

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. CC BY 4.0 (`cc-by-4.0`).
