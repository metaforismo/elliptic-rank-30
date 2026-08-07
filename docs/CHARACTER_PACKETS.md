# Galois-character packets

For a multiquadratic extension `L/K` with squareclass space `D`, the rational
Mordell--Weil space decomposes into quadratic character channels:

```text
E(L) tensor Q = direct_sum_{d in D} E^d(K) tensor Q.
```

Distinct channels are height-orthogonal.  Search records must compute every
channel, including hidden product twists, and rank the packet by total exact
character rank rather than by the number of multisections originally used to
construct it.
