# 02: Cell Growth Kinetics

## Growth phases

1. **Lag phase**: cells attach, spread and adapt after seeding or passage.
2. **Log (exponential) phase**: constant specific growth rate. Measure kinetics here.
3. **Stationary phase**: contact inhibition at confluence, or nutrient depletion.
4. **Decline**: death outpaces division when medium is exhausted.

## Exponential growth

$$\frac{dN}{dt} = \mu N \quad\Rightarrow\quad N(t) = N_0 e^{\mu t}$$

- Specific growth rate from two counts: $\mu = \dfrac{\ln(N_2/N_1)}{t_2 - t_1}$
- Doubling time: $t_d = \dfrac{\ln 2}{\mu}$
- Population doublings: $PD = \log_2(N_f / N_i)$. Sum PD across passages to
  get cumulative population doublings (CPD), the usual measure of culture age.

## Bounded growth

**Logistic**, with carrying capacity $K$ (e.g. cells at confluence):

$$\frac{dN}{dt} = rN\left(1 - \frac{N}{K}\right) \quad\Rightarrow\quad N(t) = \frac{K}{1 + \frac{K - N_0}{N_0} e^{-rt}}$$

**Gompertz**, where the growth rate itself decays (often fits spheroids and tumours):

$$N(t) = K \exp\!\left[\ln\!\left(\frac{N_0}{K}\right) e^{-rt}\right]$$

**Monod**, substrate-limited growth (from microbial kinetics; also used for
glucose- or oxygen-limited growth):

$$\mu = \mu_{max} \frac{S}{K_s + S}$$

## Worked example

5.0 × 10⁵ cells are seeded; 72 h later you count 3.2 × 10⁶ (all in log phase).

- $\mu = \ln(3.2\times10^6 / 5.0\times10^5)/72 = \ln 6.4 / 72 = 0.0258\ \text{h}^{-1}$
- $t_d = 0.693 / 0.0258 = 26.9\ \text{h}$
- $PD = \log_2 6.4 = 2.68$ doublings (and $72 / 2.68 = 26.9$ h, a consistency check)

```python
from tissue_eng import growth
mu = growth.specific_growth_rate(5e5, 3.2e6, 72)   # 0.0258 1/h
growth.doubling_time(mu)                           # 26.9 h
growth.population_doublings(5e5, 3.2e6)            # 2.68
```

## Practical points

- Fit growth curves with log-transformed counts. A straight line on a
  semi-log plot confirms exponential growth; only use those points for μ.
- Report n (biological replicates) and the spread of the counts, not just the mean.
- Seeding density changes lag time and time-to-confluence. Too sparse and
  some primary cells will not proliferate at all.
