# Oosawa Model + Dissociation — Gillespie Simulation

**Course:** Physics of Molecular Diseases — Week 1  
**Question 3:** Investigate how adding a dissociation term affects aggregation characteristics.

---

## Overview

This script extends the Oosawa model by adding a **dissociation reaction**: polymers can lose one monomer at a time. This introduces a competition between growth (association) and shrinkage (dissociation), leading to a genuine **steady state** — something the pure Oosawa model cannot produce.

The key question is: under what conditions do polymers survive, and what do they look like at steady state?

---

## The Model

### Reactions

| Reaction | Description | Rate |
|----------|-------------|------|
| **Nucleation** | $n_c$ monomers → polymer of size $n_c$ | $k_n \cdot m(t)^{n_c}$ |
| **Association** | polymer$(j)$ + monomer → polymer$(j+1)$ | $k_a \cdot m(t)$ per polymer |
| **Dissociation** | polymer$(j)$ → polymer$(j-1)$ + monomer | $k_d$ per polymer |

Note: dissociation rate does **not** depend on $m(t)$ — a monomer leaves the polymer regardless of the environment. If a polymer shrinks below $n_c$, it dissolves completely and releases all its monomers back to the free pool.

### Moment equations (ODE theory)

$$\frac{dP}{dt} = k_n \cdot m(t)^{n_c}$$

$$\frac{dM}{dt} = \left(k_a \cdot m(t) - 2k_d\right) P(t) + n_c k_n m(t)^{n_c}$$

with $m(t) = m_0 - M(t)$.

The dissociation term enters only in $\frac{dM}{dt}$, not in $\frac{dP}{dt}$ — because dissociation shrinks polymers but does not destroy them (unless they fall below $n_c$), so it has negligible effect on the polymer count for large polymers.

---

## Steady State Analysis

Setting $\frac{dM}{dt} = 0$ and neglecting the nucleation term (small at steady state):

$$\left(k_a \cdot m^* - 2k_d\right) P^* = 0$$

This gives the **critical monomer concentration**:

$$\boxed{m^* = \frac{2k_d}{k_a}}$$

### Meaningful constraints for steady state

| Condition | Outcome |
|-----------|---------|
| $m_0 > m^* = \dfrac{2k_d}{k_a}$ | Polymers grow and reach steady state at $M^* = m_0 - m^*$ |
| $m_0 < m^* = \dfrac{2k_d}{k_a}$ | Dissociation dominates — no stable polymers form |
| $k_d = 0$ | Pure Oosawa: all monomers eventually polymerise, no steady state |

The steady state polymer mass is:

$$M^* = m_0 - m^* = m_0 - \frac{2k_d}{k_a}$$

### Effect of parameters on steady state

| Parameter change | $m^*$ | $M^*$ | $f(j)$ distribution |
|-----------------|-------|-------|----------------------|
| ↑ $k_d$ | ↑ higher | ↓ lower | Narrower, smaller polymers |
| ↑ $k_a$ | ↓ lower | ↑ higher | Wider, larger polymers |
| ↑ $k_n$ | unchanged | unchanged | More polymers, same mean size |

---

## Effect on Size Distribution f(j)

Without dissociation (pure Oosawa), $f(j)$ is monotonically decreasing — dominated by the smallest polymers (sin-like shape).

Adding dissociation **truncates the tail**: very large polymers are unstable because dissociation events accumulate. At steady state, $f(j)$ becomes exponential $(\sim e^{-j})$ — still decreasing, but more sharply cut off at large $j$. The higher $k_d$, the more compressed the distribution.

A **bell-shaped (peaked) distribution** requires fragmentation, not dissociation — as discussed in the lectures.

---

## Implementation note: why t_max instead of m == 0

In the pure Oosawa model, the simulation ends when $m = 0$ (all monomers consumed). With dissociation, $m(t)$ converges to $m^* > 0$ and **never reaches zero**. Using `m == 0` as a stopping condition would cause the simulation to run indefinitely.

The fix is to stop by simulation time instead:

```python
# Pure Oosawa: stop when monomers run out
if m <= 0:
    break

# With dissociation: stop at a fixed time t_max
if t >= t_max:
    break
```

Choose `t_max` large enough that the system has clearly reached steady state.

---

## Background

Dissociation is the simplest mechanism that introduces **reversibility** into polymer growth. In the context of protein misfolding diseases, this corresponds to the dynamic exchange of monomers at fibril ends — experimentally measurable via techniques such as seeding assays and hydrogen-deuterium exchange. The critical concentration $m^*$ is a thermodynamic quantity related to the free energy difference between the monomer and polymer states, making it a key target for therapeutic intervention.

**Key references:**  
- Prof. Ala Trusina, Lecture notes from the course "Physics of Molecular Diseases", Niels Bohr Institute, 2020
- Oosawa & Asakura (1975) — original polymer nucleation-elongation model  
- Knowles et al., *Science* (2009) — nucleated polymerisation framework  
