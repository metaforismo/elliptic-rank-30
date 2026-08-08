# Exact dessin enumeration for the maximal polynomial level

**Status:** finite combinatorial computation.  The generated certificate counts
all topological dessins for the maximal nondegenerate degree-(4,7) polynomial
ansatz.  It does not yet reconstruct the Belyi maps arithmetically.

## Passport

The Mason--Stothers equality case gives a degree-14 Belyi map

\[
f=\frac{Q^2}{R},\qquad f-1=\frac{v^2L^3}{R},
\]

with cycle partitions

\[
\sigma_0:(2^7),\qquad
\sigma_1:(3^4,2),\qquad
\sigma_\infty:(11,1,1,1). \tag{1}
\]

The permutation convention in the verifier is

\[
\sigma_0\sigma_1\sigma_\infty=1.
\]

## Exhaustive algorithm

The program fixes

\[
\sigma_\infty=(0\ 1\ \cdots\ 10)(11)(12)(13).
\]

A permutation of type \((2^7)\) is a perfect matching of 14 labels.  There are
exactly

\[
13!!=135135
\]

such matchings.  For each one the program computes

\[
\sigma_1=\sigma_0\sigma_\infty^{-1},
\]

checks the target cycle partition, and checks transitivity.  Isomorphism
classes are obtained by simultaneous conjugacy under the complete centralizer
of the fixed \(\sigma_\infty\):

\[
C_{S_{14}}(\sigma_\infty)\cong C_{11}\times S_3,
\]

of order \(66\).  Lexicographic canonicalization under all 66 elements removes
all duplicate labellings.  The stabilizer in this centralizer gives the exact
automorphism-group order of each dessin.

The generated certificate records:

- the number of labelled passport triples;
- the number of transitive labelled triples;
- the number of dessin isomorphism classes;
- canonical cycle representatives;
- automorphism-group orders;
- a canonical SHA-256 digest.

## Reproduction

```bash
python3 research/degree47_belyi_passport.py \
  --output certificates/degree47_belyi_passport.json
python3 research/degree47_belyi_passport.py \
  --check certificates/degree47_belyi_passport.json
```

The computation uses only the Python standard library and enumerates the full
finite search space; there is no random sampling.

## Arithmetic next step

A topological dessin is not yet a rational polynomial decomposition.  For each
canonical triple the next stage must determine:

1. field of moduli and field of definition;
2. exact Belyi map with passport (1);
3. a normalization placing the unique order-two point over 1 at \(v=0\);
4. whether the denominator can be written
   \(v^3-Sv^2+3v+1\) over \(\mathbf Q\);
5. whether the resulting \(L,Q\) produce a non-torsion section and an
   independent rank contribution.

The two repeated-root values \(S=-3,15/4\) are outside this generic passport
and remain separate exact systems.
