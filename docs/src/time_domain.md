# Time-domain migration

[`GasChromatographySimulator.simulate`](@ref) integrates along the column. Position runs from the inlet `x = 0` to the outlet `x = L`, and the unknowns are the time `t(x)` and the temporal peak variance `τ²(x)`. That is the default, and it is the solver to use for a vacuum outlet.

The functions on this page integrate in laboratory time. The unknowns are the band position `x(t)` and the spatial variance `σ²(t)`. Use them when both end pressures stay finite and the band may stop or travel back toward `x = 0`:

- a pause, `pin = pout`, where the band stays in place
- a reversal, `pin < pout`, where the band leaves through `x = 0`
- a capillary trap that is loaded in one direction and emptied in the other

`Options` does not switch `simulate` onto these equations. Call the time-domain functions yourself.

The model is the same linear isotherm and equilibrium retention factor as the position-domain solver. It predicts where a trace-level molecule sits and when it leaves. It does not predict the sample volume that saturates a trap.

## Equations

The solute velocity is the mobile-phase velocity reduced by retention:

``u_s(x,t) = \frac{u_M(x,t)}{1+k(x,t)}, \qquad \frac{dx}{dt} = u_s``

The spatial variance obeys

``\frac{d\sigma^2}{dt} = \mathcal{D}(x,t) + 2\sigma^2\frac{\partial u_s}{\partial x}``

``\mathcal{D}`` is the Golay plate height times the solute speed, expanded so that it stays finite when the band is parked. At `u_M = 0` only longitudinal diffusion remains, ``\mathcal{D} = 2 D_M/(1+k)``. The retention factor is not clamped by `k_th` on this path: `u_s ≈ 0` is a normal state.

At an outward exit the temporal width is ``\tau = \sigma/|u_s|``. While the band is parked, `u_s = 0` and `τ` is not a meaningful width; keep `σ`.

## When to stay with the position-domain solver

A vacuum outlet, `pout = 0`, cannot reverse the flow. `pin² − pout²` stays positive for every `pin > 0`. At `x = L` the local pressure is zero, so `u_s` and `σ²` grow without a bound. The residency `r = 1/u_s` goes to zero there, which is why `dt/dx` and `dτ²/dx` stay regular. Use [`GasChromatographySimulator.solving_migration`](@ref), [`GasChromatographySimulator.solving_peakvariance`](@ref), or [`GasChromatographySimulator.simulate`](@ref) for that program.

## Initial values

[`GasChromatographySimulator.Substance`](@ref) carries two extra fields, both defaulting to `0`:

- `x₀`: position at the start, in m. `0` is the `x = 0` end.
- `σ₀`: spatial band width at the start, in m (standard deviation).

`t₀` is still the start time. `τ₀` and `σ₀` are independent. The position-domain solvers read `τ₀` and always start at `x = 0`. The time-domain solvers read `x₀` and `σ₀`. If the caller has a temporal injection width and the band is moving, convert once at the call: ``\sigma_0 = \tau_0\,|u_s(x_0, t_0)|``.

The default end time is the last absolute program time, `prog.T_itp.grid[2][end]`. Pass `t_end` to stop earlier. Program knots strictly between `t₀` and `t_end` are passed to the solver as `tstops`.

## Solvers

[`GasChromatographySimulator.solving_migration_time`](@ref) integrates `x(t)` only. In the solution, `sol.t` is time and `sol.u` is position.

[`GasChromatographySimulator.solving_bandvariance_time`](@ref) integrates `σ²(t)` along that migration. Its time span is the span of the migration solution, so the variance is not continued after the band has left the column.

[`GasChromatographySimulator.solving_odesystem_time`](@ref) integrates both as one system. `sol.u` is the pair `(x, σ²)`.

An outward exit terminates the solve: `x = L` while `u_s > 0`, or `x = 0` while `u_s < 0`. The return code is then `Terminated`, and `sol.t[end]` is the exit time. If the band is still inside at `t_end`, the return code is `Success` and `sol.u[end]` is the position (and `σ²`, for the combined system) at that time.

[`GasChromatographySimulator.peaklist`](@ref) and [`GasChromatographySimulator.sol_extraction`](@ref) read a position-domain solution, where `sol.t` is position. Do not pass a time-domain solution to them.

## Example

One isothermal column, `pin > pout`, started at the inlet with `x₀ = 0` and `σ₀ = 0`. The combined system stops at the outlet.

```@example td
using GasChromatographySimulator

col = GasChromatographySimulator.Column(1.0, 0.25e-3, 0.25e-6, "DB5", "He")
prog = GasChromatographySimulator.Program(
    [0.0, 60.0, 120.0],
    [100.0, 100.0, 100.0],
    [150000.0, 150000.0, 150000.0],
    [120000.0, 120000.0, 120000.0],
)
sub = GasChromatographySimulator.Substance(
    "test", "1-1-1", 400.0, 30.0, 100.0, 1e-3, "test", 1e-6, 0.0, 0.0,
)
opt = GasChromatographySimulator.Options()
sol = GasChromatographySimulator.solving_odesystem_time(col, prog, sub, opt)
us = GasChromatographySimulator.solute_velocity(col.L, sol.t[end], col, prog, sub, opt)
(tR = sol.t[end], x = sol.u[end][1], τ = sqrt(sol.u[end][2]) / abs(us))
```

`tR` is the exit time in seconds, `x` is the outlet, and `τ` is the temporal width there. The same column with `pout = 0` belongs on [`GasChromatographySimulator.solving_odesystem_r`](@ref) instead.
