# Segment handoff

How you pass a solute from one simulated piece of the GC path to the next depends on the **connection**, not on whether you integrate in time or in position. Choose the physical quantities first; then map them into the next solver’s state.

Programs on connected segments must share the **same absolute time axis** (temperature and pressure interpolants evaluated at the same laboratory clock).

## Sequential capillaries (junction)

Typical cases: analytical column → restrictor, column → modulator → second column, or any graph edge in [GasChromatographySystems.jl](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl).

The next inlet only sees a **time** profile of arriving mass. From the previous peak list (or exit record):

| Next initial | From previous |
|---|---|
| ``t_0`` | retention / exit time ``t_R`` (centre) |
| ``\tau_0`` | temporal width ``\tau_R`` at that exit |

Start the next segment at its own inlet (``x = 0`` in that segment’s coordinate).

**Do not** copy the previous spatial band width ``\sigma_R`` across a junction when velocity, diameter, temperature, or film change. Spatial and temporal widths are related by the local solute speed,

``\sigma = |u_s|\,\tau``,

so the same temporal peak maps to a different ``\sigma`` on each side of the junction.

### Mapping into each solver

- **Position-domain** ([`simulate`](@ref), [`solving_odesystem_r`](@ref)): set [`Substance`](@ref) fields ``t_0`` and ``\tau_0``. Integration starts at ``z = 0``.
- **Time-domain** ([`solving_odesystem_time`](@ref)): set ``t_0``, ``x_0 = 0``, and convert once at the inlet while the band is moving:

  ``\sigma_0 = \tau_0\,|u_s(0, t_0)|``.

  Keep ``\tau_0`` for bookkeeping if you like; the time-domain ODE reads ``\sigma_0``.

An open end or a focused reinjection (trap release, valve slice) is still a **junction**: estimate the temporal window of mass that leaves (e.g. from mass-fraction times ``t_{05}``–``t_{95}``, or from ``\tau_R``), then hand off ``(t_0, \tau_0)`` as above. A snapshot of trap-outlet ``\sigma`` while the band is still exiting is not the next column’s ``\sigma_0``.

## Same capillary, changing conditions

Typical cases: load → hold → reverse flush on one trap capillary; a pause (`pin = pout`) then resumed flow on the same column.

The band never leaves the tube. Carry the **spatial** state:

| Next initial | From previous end |
|---|---|
| ``x_0`` | position ``x`` |
| ``\sigma_0`` | spatial width ``\sigma`` |
| ``t_0`` | absolute time at the end of the previous phase (if the next program uses the same clock), or `0` with a phase-local program |

Truncate to ``[0, L]`` when enforcing an open-end mass balance. While parked (``u_s = 0``), ``\tau`` is not a meaningful width; keep ``\sigma``.

### Mapping into each solver

- **Time-domain**: use ``(t_0, x_0, \sigma_0)`` directly.
- **Position-domain**: only if the next phase starts at the segment inlet with forward elution; convert ``\tau_0 = \sigma_0 / |u_s|`` at that inlet. A parked mid-column plug is not a natural position-domain initial condition—prefer the time-domain solver.

## Summary

| Connection | Pass | Not this |
|---|---|---|
| Capillary → next capillary (junction) | ``t_R \to t_0``, ``\tau_R \to \tau_0`` (then ``\sigma_0 = \tau_0|u_s|`` if time-domain) | previous ``\sigma_R`` as next ``\sigma_0`` |
| Same capillary, new program phase | ``x``, ``\sigma`` (and absolute ``t``) | resetting to ``x_0 = 0`` with a copied ``\tau`` unless the band really re-enters at the inlet |

See also [Time-domain migration](@ref) for ``(x,\sigma^2)`` equations and when to stay with the position-domain solver (especially vacuum outlets).
