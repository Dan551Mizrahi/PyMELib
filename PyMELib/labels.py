from enum import IntEnum

class Label(IntEnum):
    """
    Labels for the preprocessing and enumeration of minimal dominating sets for further details see our paper:
    Batya Kenig and Dan Shlomo Mizrahi. Enumeration of minimal hitting sets parameterized by treewidth.
    In Proceedings of the 28th International Conference on Database Theory (ICDT), 2025

    In short, the labels are as follows:
    SIGMA_I (SI): Vertices in D with no private neighbor in V(G)\D and N(v) ∩ D = ∅.
    SIGMA_0 (S0): Vertices in D with a private neighbor in V(G)\D and no private neighbors in V(G)\Vi.
    SIGMA_1 (S1): Vertices in D with a private neighbor in V(G)\D .
    OMEGA_0 (W0): Vertices in V(G)\D that are private neighbors of D with a single neighbor in Vi ∩ D.
    OMEGA_1 (W1): Vertices in V(G)\D that are private neighbors of D with a single neighbor in V(G)\Vi ∩ D.
    RHO_0 (R0): Vertices in V(G)\D with |N(v) ∩ D| ≥ 2 and |N(v) ∩ Vi ∩ D| ≥ 2.
    RHO_1 (R1): Vertices in V(G)\D with |N(v) ∩ D| ≥ 2 and |N(v) ∩ Vi ∩ D| ≥ 1.
    RHO_2 (R2): Vertices in V(G)\D with |N(v) ∩ D| ≥ 2 and |N(v) ∩ Vi ∩ D| ≥ 0.
    """

    class F_sigma(IntEnum):
        SI = 0
        S0 = 1
        S1 = 2

    class F_omega(IntEnum):
        W0 = 3
        W1 = 4

    class F_rho(IntEnum):
        R0 = 6
        R1 = 7
        R2 = 8

    @property
    def in_sigma(self) -> bool:
        return self < 3

    @property
    def in_omega(self) -> bool:
        return 2 < self < 5

    @property
    def in_rho(self) -> bool:
        return 4 < self

    def same_class(self, other) -> bool:
        return self // 3 == other // 3
