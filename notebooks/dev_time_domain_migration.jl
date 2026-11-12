### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 423a3ec8-2ed1-46a8-9c2d-57b35d3db73f
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.develop(path="/Users/janleppert/Documents/GitHub/GasChromatographySimulator")
    Pkg.add([
        Pkg.PackageSpec(name="Plots", version="1"),
        Pkg.PackageSpec(name="PlutoUI", version="0.7"),
        Pkg.PackageSpec(name="SciMLBase"),# version="3"),
        Pkg.PackageSpec(name="ForwardDiff", version="1"),
    ])
    using GasChromatographySimulator, Plots, PlutoUI, ForwardDiff
    using SciMLBase
end

# ╔═╡ 31e1003e-cc5d-11f1-ac55-13ff6939a21b
begin
	#using GasChromatographySimulator
	#using PlutoUI
	#using Plots
	#import SciMLBase
	#using SciMLBase: VectorContinuousCallback, terminate!, solve
	#using ForwardDiff
end

# ╔═╡ a1e55ea0-6417-4009-931d-8ce2a2a48ecc
# ╠═╡ disabled = true
#=╠═╡
begin
	import Pkg
Pkg.develop(path="/Users/janleppert/Documents/GitHub/GasChromatographySimulator")
end
  ╠═╡ =#

# ╔═╡ 09a22481-7dba-460a-b7db-e15bd37607e8
TableOfContents(depth=3)

# ╔═╡ 56263d3f-dfdc-4863-b16a-decb6ae5c078
md"""
# Development for time-domain migration

Leave the existing solver alone. Do not change simulate, solving_migration, solving_peakvariance, or peakode. The new functions are a second entry point. peakode still returns `H*r²` while the paper says `H*r`; record that later, do not “correct” it in this branch.

1. **Add `mobile_phase_velocity` in src/Model.jl**, next to flow and mobile_phase_residency. Pressure control only. Build it from the signed normalized flow, `u_M ∝ F·T/(p·d²)`, with `F ∝ pin² − pout²`. When `pin² ≈ pout²`, return 0 instead of dividing by the pressure difference. The existing pressure profile is already valid for `pin < pout`.
Test that function before any ODE. One isothermal column, constant d and df, is enough, and the same fixture can be reused for every later test.
`pin > pout`: `u_M` matches `1/mobile_phase_residency`.
`pin < pout`: `u_M` is the negative of the swapped-end case.
`pin = pout`: `u_M` is 0.

2. **Add `solute_velocity` in the same file**. `u_s = u_M/(1+k)` using the existing `retention_factor`, and do not apply `k_th` on this path. Check the sign, and check that a large `k` makes `|u_s|` tiny but finite.

3. **Add `spatial_variance_rate` in src/Model.jl**. Evaluate the Golay terms from `plate_height` at `|u_M|`, then multiply by `|u_s|`. Do not call `plate_height` and multiply by `u_M`, because that residency-based velocity becomes negative in reverse flow. Check that the rate is never negative, and that as `u_M → 0` it tends to `2 D_M/(1+k)`.

4. **Add `solving_migration_time` in src/Solving.jl with position only**. State `x`, independent variable `t`, right-hand side `solute_velocity`. Pass `tstops` at `prog.time_steps`. A `ContinuousCallback` terminates only on an outward crossing: `x = L` while `u_s > 0`, or `x = 0` while `u_s < 0`. No exit before the end of `tspan` is a normal result. Do not copy the Pneumatics clamp that blocks `x = 0`.
Compare one forward run with `solving_migration`. Isothermal, `pin > pout`, moderate `k`. The exit time at `x = L` should match `t(L)` from the position-domain solve. This is the gate for the rest of the branch.

5. **Then add `solving_odesystem_time` with state `(x, σ²)`**. Same callbacks. `σ₀` is an argument. After test for `solving_migration_time` and `solving_migration` is green, add the reverse-flow exit, the `pin = pout` pause (`x` stays put, `σ²` does not fall), and the parked-band diffusion check.
Only after those tests, compare temporal width on a strictly forward run. Translate the published equation, `dτ²/dt = H + 2τ² u_s ∂r/∂t` (`dσ²/dt = 𝒟 + 2σ²∂u_s/∂x`), and compare it with `solving_peakvariance`. Write down the gap to `peakode`. That comparison is a note, not a change to the old solver.

6. **Add new parameters to substructures**. Add `x₀`, `σ₀` to `Substance` structure and `t_end` to `Program` structure. 

7. **Docstrings in src/Model.jl and src/Solving.jl** once the tests pass. docs/src/functions.md and docs/src/usage.md can follow. The experimental trap notebook waits until solving_odesystem_time is on this branch.
"""

# ╔═╡ 89923dd0-3ac6-4b73-ae22-4fdf16c5e9b1
md"""
## ✅ 1. `mobile_phase_velocity`
"""

# ╔═╡ c6ddea2e-e836-4787-bbe2-6b6892f0366d
function mobile_phase_velocity(x, t, T_itp, Fpin_itp, pout_itp, L, d::Number, gas; ng=false, vis="Blumberg", control="Pressure")
	if control == "Pressure"
		pin_itp = Fpin_itp
		pp = GasChromatographySimulator.pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control="Pressure")
		κL = GasChromatographySimulator.flow_restriction(L, t, T_itp, d, gas; ng=ng, vis=vis)
		uM = 1/64*T_itp(x, t)/(κL*d^2)*(pin_itp(t)^2-pout_itp(t)^2)/pp
	elseif control == "Flow"
		F_itp = Fpin_itp
		pp = GasChromatographySimulator.pressure(x, t, T_itp, F_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control="Flow")
		uM = 4/π * GasChromatographySimulator.pn/GasChromatographySimulator.Tn * F_itp(t)/(d^2) * T_itp(x,t)/pp
	end
return uM # in m/s
end

# ╔═╡ d508f096-23d8-4275-a7dc-75756b0ebaec
function mobile_phase_velocity(x, t, par::GasChromatographySimulator.Parameters)
    uM = mobile_phase_velocity(x, t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
    return uM
end

# ╔═╡ 8306a4d2-d3ec-44fc-b384-0dfec29502bb
md"""
### ✅ Test
"""

# ╔═╡ cd5b925a-15d9-4470-a9c1-5d6bb27979c3
col = GasChromatographySimulator.Column(1.0, 0.25e-3, 0.25e-6, "DB5", "He")

# ╔═╡ f98a6634-119b-4c57-ba52-336f08c3f354
prog1 = GasChromatographySimulator.Program([0.0, 60.0, 120.0], [100.0, 100.0, 100.0], [150000.0, 150000.0, 150000.0], [120000.0, 120000.0, 120000.0])

# ╔═╡ 67b656fa-deaa-454a-9954-2c4ca61884d1
prog1.T_itp.grid[2]

# ╔═╡ 04ff0803-eb2e-49fb-977d-742093f201f3
sub = GasChromatographySimulator.Substance("test1", "1-1-1", 400.0, 30.0, 100.0, 1e-3, "test", 1e-6, 0.0, 0.0)

# ╔═╡ ebb32727-9603-4fe9-9a06-42ba80a16fff
opt = GasChromatographySimulator.Options()

# ╔═╡ b297bbbf-eb7d-4353-9a5e-e59f30af34d4
par1 = GasChromatographySimulator.Parameters(col, prog1, [sub], opt)

# ╔═╡ 9cec64ad-3645-4447-bd93-391890bc0422
# pressure reversed
prog2 = GasChromatographySimulator.Program([0.0, 60.0, 120.0], [100.0, 100.0, 100.0], [120000.0, 120000.0, 120000.0], [150000.0, 150000.0, 150000.0])

# ╔═╡ 195ed35a-61c9-43ea-807e-a8c7f187f8c0
par2 = GasChromatographySimulator.Parameters(col, prog2, [sub], opt)

# ╔═╡ 7388cf27-8603-454c-bd9e-92469a7fda02
uM1(x,t) = mobile_phase_velocity(x, t, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.col.L, par1.col.d, par1.col.gas; ng=par1.opt.ng, vis=par1.opt.vis, control=par1.opt.control)

# ╔═╡ 049a20a9-5b86-467a-a752-c022af50ec6c
uM1(0.0, 0.0), uM1(1.0, 0.0)

# ╔═╡ 0b10c438-0104-4ba6-9491-9610aa8e5181
uM2(x,t) = mobile_phase_velocity(x, t, par2)

# ╔═╡ ddd701f6-b8e9-4cf6-b986-f4f631d456b6
uM2(0.0, 0.0), uM2(1.0, 0.0)

# ╔═╡ 7153bf12-e64b-4679-8538-1ec44f151fa4
uM1(0.0, 0.0) == - uM2(1.0, 0.0)

# ╔═╡ 85edb9e5-ee3a-4e36-acb9-2dd036584730
uM1(1.0, 0.0) == - uM2(0.0, 0.0)

# ╔═╡ 74f14b83-3e25-4813-b943-bbf724814f0e
rM1(x, t) = GasChromatographySimulator.mobile_phase_residency(x, t, par1)

# ╔═╡ 8dbae989-f522-456e-a051-6bfe34eebabe
rM1(0.0, 0.0), rM1(1.0, 0.0)

# ╔═╡ 38316056-091c-4097-a551-7439f9b4c995
uM1(0.0, 0.0) == 1/rM1(0.0, 0.0)

# ╔═╡ d30479bd-8fc2-4990-b1ee-cc5d506e8d7c
uM1(1.0, 0.0) == 1/rM1(1.0, 0.0)

# ╔═╡ 27e86e63-ed5c-44ce-bad8-74ae2f55c21c
uM1(1.0, 0.0) - 1/rM1(1.0, 0.0)

# ╔═╡ de03478f-7af7-40db-8e57-e09552db1667
uM1(1.0, 0.0) ≈ 1/rM1(1.0, 0.0)

# ╔═╡ aa5d9359-8771-4c0a-ac98-4dc5e0b34fdc
rM2(x, t) = GasChromatographySimulator.mobile_phase_residency(x, t, par2)

# ╔═╡ 78cae962-2eb3-4bb0-a6a7-62457f16f85d
rM2(0.0, 0.0), rM2(1.0, 0.0)

# ╔═╡ be4cd698-2ba3-4666-bbbd-e91cd9864e47
uM2(0.0, 0.0) == 1/rM2(0.0, 0.0)

# ╔═╡ 0d5946bd-2fad-4bcd-8bb6-264b4d895a1d
uM2(0.0, 0.0) - 1/rM2(0.0, 0.0)

# ╔═╡ 0655e3d3-f4cb-4224-8c05-64d7cdb9e3e3
uM2(0.0, 0.0) ≈ 1/rM2(0.0, 0.0)

# ╔═╡ 07875f26-915e-41b5-8aff-589d891c4f99
uM2(1.0, 0.0) == 1/rM2(1.0, 0.0)

# ╔═╡ 52ec2ab2-3f37-48d9-be9e-73289ed86ff3
# same pressure
prog3 = GasChromatographySimulator.Program([0.0, 60.0, 120.0], [100.0, 100.0, 100.0], [150000.0, 150000.0, 150000.0], [150000.0, 150000.0, 150000.0])

# ╔═╡ 9fd00d96-0589-403f-bbf3-2bc6ed7b1e24
opt3 = GasChromatographySimulator.Options(k_th=Inf)

# ╔═╡ 5d154824-9280-4d7f-9758-58d96fc4136b
par3 = GasChromatographySimulator.Parameters(col, prog3, [sub], opt3)

# ╔═╡ 3262fb55-48b4-407b-a844-df8248afa789
uM3(x,t) = mobile_phase_velocity(x, t, par3)

# ╔═╡ 7beafde2-7752-4dff-a177-4c2e0cad8a1f
uM3(0.0, 0.0), uM3(1.0, 0.0)

# ╔═╡ 768295ef-1658-4ab6-ba85-59945cb8a8dd
rM3(x, t) = GasChromatographySimulator.mobile_phase_residency(x, t, par3)

# ╔═╡ 8a7e4645-9504-415a-a345-cce6ef5c1808
rM3(0.0, 0.0), rM3(1.0, 0.0)

# ╔═╡ f1bd774e-137f-4895-8517-e4f87252649c
md"""
#### Defined flow
"""

# ╔═╡ 15f972bb-58d4-49b1-b321-a2f82f322fd4
progF = GasChromatographySimulator.Program([0.0, 60.0, 120.0], [100.0, 100.0, 100.0], [1/60e6, 1/60e6, 1/60e6], [120000.0, 120000.0, 120000.0])

# ╔═╡ 1a8fcf4d-9cf5-41c2-97dc-78c2a66f126a
optF = GasChromatographySimulator.Options(control="Flow")

# ╔═╡ a7bac514-b370-4ed2-8b0e-0cbfd9b3b0bc
parF = GasChromatographySimulator.Parameters(col, progF, [sub], optF)

# ╔═╡ 7b461249-425f-4bf9-ab76-d419ce61ccab
uMF(x,t) = mobile_phase_velocity(x, t, parF)

# ╔═╡ a0964a38-0124-47f3-9f12-4f25a25007f6
uMF(0.0,0.0), uMF(1.0,0.0)

# ╔═╡ dbcfa37e-02f8-4b1f-9871-6936e53052be
progFrev = GasChromatographySimulator.Program([0.0, 60.0, 120.0], [100.0, 100.0, 100.0], [-1/60e6, -1/60e6, -1/60e6], [120000.0, 120000.0, 120000.0])

# ╔═╡ b8282421-962d-4ee2-b65e-8b6f260bde1e
md"""
## ✅ 2. `solute_velocity`
"""

# ╔═╡ ff76d100-9c5b-4068-93e3-62790f45549a
function solute_velocity(x, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=false, vis="Blumberg", control="Pressure")
	uM(x, t) = mobile_phase_velocity(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)
	k(x, t) = GasChromatographySimulator.retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀; k_th=Inf)	
	return uS = uM(x, t)/(1+k(x,t))
end

# ╔═╡ e839ccda-04cc-4f2b-9a60-248a230b17d2
function solute_velocity(x, t, col::GasChromatographySimulator.Column, prog::GasChromatographySimulator.Program, sub::GasChromatographySimulator.Substance, opt::GasChromatographySimulator.Options)
	uS = solute_velocity(x, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=opt.ng, vis=opt.vis, control=opt.control)
	return uS
end

# ╔═╡ 62533e44-c765-4a54-932d-0a6d6b2024d6
md"""
### ✅ Test
"""

# ╔═╡ d6c17996-1ead-4281-bc7e-7629b0afb745
uS1(x,t) = solute_velocity(x, t, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.col.L, par1.col.d, par1.col.df, par1.col.gas, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀; ng=par1.opt.ng, vis=par1.opt.vis, control=par1.opt.control)

# ╔═╡ fb02ba01-0e4f-4e1a-8bc1-5e805a0e1f40
uS2(x,t) = solute_velocity(x, t, par2.col, par2.prog, par2.sub[1], par2.opt)

# ╔═╡ 81a48945-ea78-4ac9-b128-c8426cd7c071
uS1(0.0, 0.0), uS1(1.0, 0.0)

# ╔═╡ 5ad305c4-dfd5-40fd-b422-cda482ab7acf
uS2(0.0, 0.0), uS2(1.0, 0.0)

# ╔═╡ dcbb2198-66eb-48d1-a4a1-54375de2bf98
uS1(0.0, 0.0) == - uS2(1.0, 0.0)

# ╔═╡ 0654d56f-2458-450f-995d-257c86abe2a3
uS1(1.0, 0.0) == - uS2(0.0, 0.0)

# ╔═╡ c08c8663-e257-4228-b358-c5f40320efb6
k(x, t) = GasChromatographySimulator.retention_factor(x, t, par1.prog.T_itp, par1.col.d, par1.col.df, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀; k_th=Inf)

# ╔═╡ fb9f99af-183b-45db-8868-5b2a199f8916
k(0.0, 0.0)

# ╔═╡ 4844421b-c7a7-4029-8dd5-dba9ab829b88
k0(x, t) = GasChromatographySimulator.retention_factor(x, t, par1.prog.T_itp, par1.col.d, par1.col.df, 0.0, 0.0, 0.0, par1.sub[1].φ₀; k_th=Inf)

# ╔═╡ 570e1699-86f9-43d9-812a-f78bf1402b54
k0(0.0, 0.0)

# ╔═╡ a43dc802-b93a-4a17-9135-67cb5a800297
khigh(x,t) = GasChromatographySimulator.retention_factor(x, t, par1.prog.T_itp, par1.col.d, par1.col.df, 1000.0, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀; k_th=Inf)

# ╔═╡ 194dc600-bb61-4040-8ae0-78c2f92ff09b
khigh(0.0, 0.0)

# ╔═╡ 3cd03859-8aaf-4ad1-a1f1-84534d42a385
uS0(x,t) = solute_velocity(x, t, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.col.L, par1.col.d, par1.col.df, par1.col.gas, 0.0, 0.0, 0.0, par1.sub[1].φ₀; ng=par1.opt.ng, vis=par1.opt.vis, control=par1.opt.control)

# ╔═╡ b7fbab24-fb3d-42fd-9c53-9e1d510ba69c
uS0(0.0, 0.0)

# ╔═╡ 2a4c8410-b5ec-4a80-95a6-22e0d264aa99
uS0(0.0, 0.0) == uM1(0.0, 0.0)

# ╔═╡ 04199e01-eac5-4b01-8a4a-504fe6b415ad
uShigh(x,t) = solute_velocity(x, t, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.col.L, par1.col.d, par1.col.df, par1.col.gas, 1000.0, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀; ng=par1.opt.ng, vis=par1.opt.vis, control=par1.opt.control)

# ╔═╡ 144dc9b7-d0ab-4f80-846c-165040f299f3
uShigh(0.0,0.0)

# ╔═╡ 0236dc1c-e1ad-480c-8e69-c087bd1085cf
uS3(x,t) = solute_velocity(x, t, par3.col, par3.prog, par3.sub[1], par3.opt)

# ╔═╡ 85f94ab5-2135-4641-ab46-b007214b91ee
uS3(0.0,0.0), uS3(1.0,0.0)

# ╔═╡ d78e12a1-b46d-4aed-80fb-9c0c4deb3983
md"""
## ✅ 3. `temporal_dispersion_rate`
"""

# ╔═╡ e51628de-5a58-4ec2-8d37-adf435b276dd
function temporal_dispersion_rate(x, t, T_itp, Fpin_itp, pout_itp, L, d::Number, df::Number, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=false, vis="Blumberg", control="Pressure")
	# temporal dispersion rate 𝒟 = spatial dispersion rate ℋ * solute velocity
    # 𝒟 is used for the time-domain solver instead of ℋ in the space-domain
	id = d# - 2.0*df
    uM = mobile_phase_velocity(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)
    μ = 1/(1 + GasChromatographySimulator.retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀, k_th=Inf))
    DM = GasChromatographySimulator.diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=ng, vis=vis, control=control)
    DS = DM/10000
    𝒟1 = 2*DM*μ
    𝒟2 = id^2/96*(6*μ^3-16*μ^2+11*μ)*uM^2/DM
    𝒟3 = 2/3*df^2*μ^2*(1-μ)*uM^2/DS
    𝒟 = 𝒟1 + 𝒟2 + 𝒟3
    return 𝒟
end

# ╔═╡ 6e78cdef-7f2f-4ed4-abce-0964905ade80
function temporal_dispersion_rate(x, t, col::GasChromatographySimulator.Column, prog::GasChromatographySimulator.Program, sub::GasChromatographySimulator.Substance, opt::GasChromatographySimulator.Options)
    𝒟 = temporal_dispersion_rate(x, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
    return 𝒟
end

# ╔═╡ 04f05dfc-4ea2-4ffe-a22c-6f06aa81e367
# modification to use the absolute value of uM
function plate_height(x, t, T_itp, Fpin_itp, pout_itp, L, d::Number, df::Number, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    id = d# - 2.0*df
    uM = 1/GasChromatographySimulator.mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)
    μ = 1/(1 + GasChromatographySimulator.retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀, k_th=k_th))
    DM = GasChromatographySimulator.diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=ng, vis=vis, control=control)
    DS = DM/10000
    ℋ1 = 2*DM/abs(uM)
    ℋ2 = id^2/96*(6*μ^2-16*μ+11)*abs(uM)/DM
    ℋ3 = 2/3*df^2*μ*(1-μ)*abs(uM)/DS
    ℋ = ℋ1 + ℋ2 + ℋ3
    return ℋ
end

# ╔═╡ 302691e9-01d9-44cc-9da6-a757fbe13012
md"""
### ✅ Test
"""

# ╔═╡ 377cc5c2-280f-4d0e-b8d5-f467a34edc99
ℋ1(x,t) = GasChromatographySimulator.plate_height(x, t, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.col.L, par1.col.d, par1.col.df, par1.col.gas, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀, par1.sub[1].Cag; ng=par1.opt.ng, vis=par1.opt.vis, control=par1.opt.control, k_th=Inf)

# ╔═╡ 00d98412-72e6-4704-9ca0-c22c1a7937bb
ℋ1(0.0,0.0), ℋ1(1.0,0.0)

# ╔═╡ 404ef427-7785-4649-bc1b-143c432c27de
ℋ2(x,t) = GasChromatographySimulator.plate_height(x, t, par2.prog.T_itp, par2.prog.Fpin_itp, par2.prog.pout_itp, par2.col.L, par2.col.d, par2.col.df, par2.col.gas, par2.sub[1].Tchar, par2.sub[1].θchar, par2.sub[1].ΔCp, par2.sub[1].φ₀, par2.sub[1].Cag; ng=par2.opt.ng, vis=par2.opt.vis, control=par1.opt.control, k_th=Inf)

# ╔═╡ 28421f8d-a107-42ac-b222-da21d5fec223
ℋ2(0.0, 0.0), ℋ2(1.0, 0.0)

# ╔═╡ b48d5a49-ae01-41da-8983-ba522ea7b811
ℋ2_(x,t) = plate_height(x, t, par2.prog.T_itp, par2.prog.Fpin_itp, par2.prog.pout_itp, par2.col.L, par2.col.d, par2.col.df, par2.col.gas, par2.sub[1].Tchar, par2.sub[1].θchar, par2.sub[1].ΔCp, par2.sub[1].φ₀, par2.sub[1].Cag; ng=par2.opt.ng, vis=par2.opt.vis, control=par1.opt.control, k_th=Inf)

# ╔═╡ 52f0f5c6-a3f2-4c71-8f8b-8fc92f7f9694
ℋ2_(0.0, 0.0), ℋ2_(1.0, 0.0)

# ╔═╡ 0cf79001-43d7-432f-ba24-cbfb6ba3b93e
𝒟1(x, t)  = temporal_dispersion_rate(x, t, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.col.L, par1.col.d, par1.col.df, par1.col.gas, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀, par1.sub[1].Cag; ng=par1.opt.ng, vis=par1.opt.vis, control=par1.opt.control, k_th=Inf)

# ╔═╡ 9480b569-3a5b-4f04-9acf-b91235120fc7
𝒟1(0.0, 0.0), 𝒟1(1.0, 0.0)

# ╔═╡ d27b6674-6857-441f-ab7a-0ea13ae83690
𝒟2(x, t)  = temporal_dispersion_rate(x, t, par2.prog.T_itp, par2.prog.Fpin_itp, par2.prog.pout_itp, par2.col.L, par2.col.d, par2.col.df, par2.col.gas, par2.sub[1].Tchar, par2.sub[1].θchar, par2.sub[1].ΔCp, par2.sub[1].φ₀, par2.sub[1].Cag; ng=par2.opt.ng, vis=par2.opt.vis, control=par2.opt.control, k_th=Inf)

# ╔═╡ 5f134571-2c33-4a2c-99e7-06050544c6ee
𝒟2(0.0, 0.0), 𝒟2(1.0, 0.0)

# ╔═╡ 9b772052-82ad-4fcc-b35e-7c9da66ce6cc
ℋ2_(0.0, 0.0)*abs(uS2(0.0, 0.0)) == 𝒟2(0.0, 0.0)

# ╔═╡ 8ca13a43-0b0a-43a5-a9d6-4e61a54299ab
ℋ2_(0.0, 0.0)*abs(uS2(0.0, 0.0)) - 𝒟2(0.0, 0.0)

# ╔═╡ a6d98683-98e5-4074-8ca6-c8eb89fd17f3
ℋ2_(0.0, 0.0)*abs(uS2(0.0, 0.0)) ≈ 𝒟2(0.0, 0.0)

# ╔═╡ d24a67e8-e01b-4b01-aea8-2aed0f5cbf13
𝒟1(0.0, 0.0) == 𝒟2(1.0, 0.0)

# ╔═╡ 2cf1c3f4-b093-4cbb-a2a4-55d429ec841e
𝒟3(x, t)  = temporal_dispersion_rate(x, t, par3.col, par3.prog, par3.sub[1], par3.opt)

# ╔═╡ a95a9a9f-e1a2-4c2a-aebc-ff100ca10c9e
𝒟3(0.0, 0.0), 𝒟3(1.0, 0.0)

# ╔═╡ 25a4a2c3-ed68-4d4a-9432-298d3a5fe27e
ℋ3_(x,t) = plate_height(x, t, par3.prog.T_itp, par3.prog.Fpin_itp, par3.prog.pout_itp, par3.col.L, par3.col.d, par3.col.df, par3.col.gas, par3.sub[1].Tchar, par3.sub[1].θchar, par3.sub[1].ΔCp, par3.sub[1].φ₀, par3.sub[1].Cag; ng=par3.opt.ng, vis=par3.opt.vis, control=par3.opt.control, k_th=par3.opt.k_th)

# ╔═╡ a8471106-8e88-4abc-8841-2ab791a72f77
ℋ3_(0.0, 0.0), ℋ3_(1.0, 0.0)

# ╔═╡ 1170bcc9-ac3d-44af-9c67-cbd9b85a7e3e
md"""
## ✅ 4. `solving_migration_time`
"""

# ╔═╡ 3fc58fe7-8387-4952-ab50-d49640bc1222
function solving_migration_time(col::GasChromatographySimulator.Column, prog::GasChromatographySimulator.Program, sub::GasChromatographySimulator.Substance, opt::GasChromatographySimulator.Options, x₀; t_end=prog.time_steps, kwargs...) # x₀ later into Substance
    solution_zt = solving_migration_time(col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, x₀, sub.t₀, opt; t_end=t_end, kwargs...)
    return solution_zt
end

# ╔═╡ f18f1d88-dbbb-4528-a803-607d60adf84b
# from GasChromatographySimulator.jl
function solving_migration(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, t₀, opt; kwargs...)
	# this version should be autodiffable for quantities in `p`
    p = (T_itp, Fpin_itp, pout_itp, L, d, df, Tchar, θchar, ΔCp, φ₀, gas, opt) # should be autodiffable for all these parameters
	f_tz(t,p,z) = residency(z, t, p[1], p[2], p[3], p[4], p[5], p[6], p[11], p[7], p[8], p[9], p[10]; ng=p[12].ng, vis=p[12].vis, control=p[12].control, k_th=p[12].k_th)
	t₀ = t₀
	zspan = (0.0, L)
	prob_tz = ODEProblem(f_tz, t₀, zspan, p)
	solution_tz = solve(prob_tz, alg=opt.alg, abstol=opt.abstol, reltol=opt.reltol; kwargs...)
    if SciMLBase.successful_retcode(solution_tz) == false
        @info "ODE solver failed, trying again with set dt=L/1000000."
        solution_tz = solve(prob_tz, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol; kwargs..., dt=L/1000000)
    end
	return solution_tz
end

# ╔═╡ a4d85b2e-197f-41f4-a602-e0bb83204cdf
# from /GasChromatographyPneumatics.jl/notebooks/DynamicPressureSystem_coupled_StevesSys.jl
"""
	    transit_time_segment(u_itp, L, t_inj, t_max; x0=0.0, t_start=0.0, mp, alg=Tsit5(),
	                          reltol=1e-6, abstol=1e-8)

	Lagrangian transit time [s] along one column: integrate ``dx/dt = u(x,t)`` from `x0`
	until `x = L` (``ContinuousCallback`` + `terminate!`).

	`u_itp(t, x)` — bilinear interpolant from `linear_interpolation((t_, x_), fields.uM)`.
	Pass `t_start` and `mp` so times beyond the saved PDE grid use periodic extrapolation.
	Inlet at `x = 0` is impermeable (no backflow out of the column).

	Returns `(; transit, t_arrival, sol)`; `transit = NaN` if `L` is not reached before
	the integration stop time.
	"""
	function transit_time_segment(u_itp, L, t_inj, t_max; x0=0.0, t_start=0.0, mp=nothing,
	                              alg=Tsit5(), reltol=1e-6, abstol=1e-8, save_everystep=false)
	    t_grid = u_itp.grid[1]
	    mp === nothing && (mp = t_grid[end] - t_grid[1])
	    t_stop = _transit_t_stop(t_inj, t_max, t_grid; mp)
	    x_lo, x_hi = extrema(u_itp.grid[2])
	    vel_time = _velocity_time_closure(t_grid, t_start, mp)

	    function f!(du, u, p, t)
	        x = clamp(u[1], x_lo, x_hi)
	        du[1] = u_itp(vel_time(t), x)
	        x <= x_lo && du[1] < 0 && (du[1] = 0.0)
	    end

	    condition(u, t, integrator) = u[1] - L
	    affect!(integrator) = SciMLBase.terminate!(integrator)
	    cb = ContinuousCallback(condition, affect!)

	    prob = ODEProblem(f!, [x0], (t_inj, t_stop))
	    sol = solve(prob, alg; callback=cb, reltol=reltol, abstol=abstol, save_everystep=save_everystep)

	    arrived = sol.retcode == SciMLBase.ReturnCode.Terminated || sol.retcode == :Terminated
	    if !arrived
	        return (; transit=NaN, t_arrival=NaN, sol)
	    end
	    t_arr = sol.t[end]
	    return (; transit=t_arr - t_inj, t_arrival=t_arr, sol)
	end

# ╔═╡ c08f9403-61bb-4478-81ec-8cd6d73cf16a
md"""
### ✅ Test
"""

# ╔═╡ 00bd501f-ca26-47ae-a39c-b48b12121b99
# 'classic' solve_migration (only for normal flow)
sol_classic = GasChromatographySimulator.solving_migration(par1.col.L, par1.col.d, par1.col.df, par1.col.gas, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀, par1.sub[1].t₀, par1.opt)

# ╔═╡ 5e7746c3-c15e-427a-a73a-fbc0120e51a2
sol_classic.u[end]

# ╔═╡ a1f284b5-5d46-44aa-b954-c768040076b6
md"""
## 5. ✅ `solving_odesystem_time`
"""

# ╔═╡ 3579dae0-c0ed-4775-a134-acde3090c50e
md"""
### ✅ `solving_bandvariance_time`
Equivalent to function `solving_peakvariance`.
"""

# ╔═╡ 71fd80e1-6784-433a-949e-a507945efdc6
function solving_bandvariance_time(solution_zt, col, prog, sub, opt; kwargs...)
    solution_σ²z = solving_bandvariance_time(solution_zt, col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag, σ₀, opt; kwargs...) # σ₀ later in Substance
    return solution_σ²z
end

# ╔═╡ 1e0aa0b0-9f68-4458-9644-1403aa036a15
function bandode(z, t, σ², L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, opt)
	𝒟 = temporal_dispersion_rate(z, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=opt.ng, vis=opt.vis, control=opt.control)

	uS(x) = solute_velocity(x, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=opt.ng, vis=opt.vis, control=opt.control)
	∂uS∂x = ForwardDiff.derivative(uS, z)

	return 𝒟 + 2*σ²*∂uS∂x
end

# ╔═╡ abd96ffd-faf4-4d34-9790-42014aefddd8
function bandode(z, t, σ², col, prog, sub, opt)
    return bandode(z, t, σ², col.L, col.d, col.df, col.gas, prog.T_itp, prog.Fpin_itp, prog.pout_itp, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag, opt)
end

# ╔═╡ 9f41ed53-5a28-4862-b172-fdd56da6612b
md"""
### ✅ Test
"""

# ╔═╡ 54eeb56a-b6f8-4bcd-a477-46cd61a0ec0b
md"""
### ✅ `solving_odesystem_time`
"""

# ╔═╡ f5db1295-221b-4cf2-8d8d-f076610c22c3
function solving_odesystem_time(
	col::GasChromatographySimulator.Column, 
	prog::GasChromatographySimulator.Program, 
	sub::GasChromatographySimulator.Substance, 
	x₀, σ₀, # latter add to Substance
	opt::GasChromatographySimulator.Options; 
	t_end=sum(prog.T_itp.grid[2]), kwargs...
)
	solution = solving_odesystem_time(
		col.L, col.d, col.df, col.gas, 
		prog.T_itp, prog.Fpin_itp, prog.pout_itp, 
		sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag, 
		sub.t₀, x₀, σ₀, 
		opt; 
		t_end=t_end, kwargs...
	)
	return solution
end

# ╔═╡ 5e7a2e30-dc78-4a8e-9284-e74e3414e898


# ╔═╡ efbf0fd3-8fb5-4292-851e-14d68472ef75
function odesystem_time!(du, u, p, t)
    L, d, df, gas = p[1], p[2], p[3], p[4]
    T_itp, Fpin_itp, pout_itp = p[5], p[6], p[7]
    Tchar, θchar, ΔCp, φ₀, Cag = p[8], p[9], p[10], p[11], p[12]
    opt = p[13]

    x = u[1]
    σ² = u[2]
    
    # Compute residency once and reuse it (performance optimization)
    du[1] = solute_velocity(x, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=opt.ng, vis=opt.vis, control=opt.control)
    du[2] = bandode(x, t, σ², L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, opt)
end

# ╔═╡ c14f63c8-3708-4d83-bbcf-596b2d88057d
md"""
### ✅ Test
"""

# ╔═╡ 8617ab5b-8657-484e-9e10-41055db08033
# 'classic' solve_migration (only for normal flow)
sol_sys_classic = GasChromatographySimulator.solving_odesystem_r(par1.col, par1.prog, par1.sub[1], par1.opt)

# ╔═╡ 7be3558a-deed-4f05-82a8-b50906b73712
sol_sys_classic.u[end][1] # retention time

# ╔═╡ 76c707f2-7b44-4058-8a3a-3af98dee3bb2
sqrt(sol_sys_classic.u[end][2]) # peak width in s

# ╔═╡ ec119d1b-0bca-4ae6-aa97-5839ae5f07dd
# near vacuum outlet
prog_vac = GasChromatographySimulator.Program([0.0, 60.0, 120.0], [100.0, 100.0, 100.0], [150000.0, 150000.0, 150000.0], [0.1, 0.1, 0.1])

# ╔═╡ 95ccf180-c5ed-4ae8-a4aa-71e7f79a7f60
par_vac = GasChromatographySimulator.Parameters(col, prog_vac, [sub], opt)

# ╔═╡ 49abc3ae-ff8a-4869-8f9d-50c3559c40cc
uS_vac(x,t) = solute_velocity(x, t, par_vac.col, 
	par_vac.prog,
	par_vac.sub[1],
	par_vac.opt)

# ╔═╡ 0a3ed9d8-40d8-47c0-87ac-84e4be66be21
uS_vac(0.0,0.0)

# ╔═╡ c214a16a-4055-4bca-91ca-9098477ea847
# 
sol_sys_classic_vac = GasChromatographySimulator.solving_odesystem_r(par_vac.col, par_vac.prog, par_vac.sub[1], par_vac.opt)

# ╔═╡ 494b6f1a-13a2-420a-958b-a120237e6092
sqrt(sol_sys_classic_vac.u[end][2])

# ╔═╡ 7cdaa4a9-eed0-4fc6-959d-755e0602fa60
md"""
## 6. `Substance` structure
"""

# ╔═╡ 64b9dd0c-7721-4df1-8a92-6b3a8541887b
struct Substance
    name::String        # name of solute
    CAS::Union{Missing,String}         # CAS registry number
    Tchar               # characteristic temperature in K
    θchar               # characteristic thermal constant in °C
    ΔCp                 # 3rd parameter
    φ₀                  # dimless film thickness for which Tchar, θchar and ΔCp were estimated
    ann::String         # annotations, e.g. the source of the data from which Tchar, θchar and ΔCp were estimated
    Cag                 # diffusitivity constant of analyt in a gas, calculate from structure (or from measurements)
    t₀                  # initial time in s  	
    τ₀                  # initial peak width in s
    x₀                  # initial start position on column in m (needs time-domain solver)
    σ₀                  # initial band width in m (needs time-domain solver)
    Substance(name, CAS, Tchar, θchar, ΔCp, φ₀, ann, Cag, t₀, τ₀, x₀, σ₀) = begin
        for (val, n) in ((Tchar, "Tchar"), (θchar, "θchar"), (ΔCp, "ΔCp"), (φ₀, "φ₀"), (Cag, "Cag"), (t₀, "t₀"), (τ₀, "τ₀"), (x₀, "x₀"), (σ₀, "σ₀"))
            isfinite(GasChromatographySimulator._value_real(val)) || throw(ArgumentError("`$n` must be finite, got `$val`."))
        end
        new(name, CAS, Tchar, θchar, ΔCp, φ₀, ann, Cag, t₀, τ₀, x₀, σ₀)
    end
end

# ╔═╡ e6e023a2-9d18-411b-8fb2-a3d682d56675
sub_new = Substance("test1", "1-1-1", 400.0, 30.0, 100.0, 1e-3, "test", 1e-6, 1.0, 0.0, 0.5, 0.1)

# ╔═╡ 3832d314-4e48-4be7-b00b-9ad588c12c54
function solving_odesystem_time(
	col::GasChromatographySimulator.Column, 
	prog::GasChromatographySimulator.Program, 
	sub::Substance, 
	opt::GasChromatographySimulator.Options; 
	t_end=sum(prog.T_itp.grid[2]), kwargs...
)
	solution = solving_odesystem_time(
		col.L, col.d, col.df, col.gas, 
		prog.T_itp, prog.Fpin_itp, prog.pout_itp, 
		sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag, 
		sub.t₀, sub.x₀, sub.σ₀, 
		opt; 
		t_end=t_end, kwargs...
	)
	return solution
end

# ╔═╡ b2905bc6-00f7-4bd6-8e44-4896ea14142f
sub_new

# ╔═╡ 100ef7fc-3162-4872-af25-fa132ae94ebc
md"""
## Helper functions
"""

# ╔═╡ 562a83c8-2373-4bce-ac89-189174573d95
function time_domain_exit_callback(L, t₀, p; vel = nothing)
    # p = (
	#	L, d, df, gas, # 1-4
	#	T_itp, Fpin_itp, pout_itp, # 5-7 
	#	Tchar, θchar, ΔCp, φ₀, Cag, # 8-12
	#	opt # 13
	#)
    # If vel is nothing, build it from that p-layout.
    uS = if vel === nothing
        (x, t) -> solute_velocity(
            x, t, p[5], p[6], p[7], p[1], p[2], p[3], p[4], p[8], p[9], p[10], p[11]; ng=p[13].ng, vis=p[13].vis, control=p[13].control,
        )
    else
        vel
    end
    function condition(out, u, t, integrator)
        x = u isa Number ? u : u[1]   # scalar x or (x, σ²)
        out[1] = x - L
        out[2] = x
    end
    function affect!(integrator, idx)
        x = integrator.u isa Number ? integrator.u : integrator.u[1]
        t = integrator.t
        us = uS(x, t)
        outward = (idx == 1 && us > 0) || (idx == 2 && us < 0 && t > t₀)
        outward && SciMLBase.terminate!(integrator)
    end
    return VectorContinuousCallback(condition, affect!, 2)
end

# ╔═╡ a75817af-3f1a-4e55-bd1c-f84a3e8a26f6
function program_knots(T_itp, t₀, t_end)
    knots = cumsum(T_itp.grid[2])
    return knots[t₀ .< knots .< t_end]
end

# ╔═╡ a3e1934a-8b28-404b-9a4f-c910bddce852
function solving_migration_time(L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, x₀, t₀, opt; t_end=sum(T_itp.grid[2]), kwargs...)

    p = (
		L, d, df, gas, # 1-4
		T_itp, Fpin_itp, pout_itp, # 5-7 
		Tchar, θchar, ΔCp, φ₀, NaN, # 8-12
		opt # 13
	)
	
	f_zt(z,p,t) = solute_velocity(z, t, p[5], p[6], p[7], p[1], p[2], p[3], p[4], p[8], p[9], p[10], p[11]; ng=p[13].ng, vis=p[13].vis, control=p[13].control)
	# initial value
    u0 = x₀
    # the time_steps used for the temperature interpolation
    time_steps = T_itp.grid[2]
	tspan = (t₀, t_end)
    knots = program_knots(T_itp, t₀, t_end)

    # callback definition for stop at either column inlet or outlet depending on flow direction
    cb = time_domain_exit_callback(L, t₀, p; vel = nothing)

	prob_zt = ODEProblem(f_zt, u0, tspan, p)
	solution_zt = solve(prob_zt; alg=opt.alg, abstol=opt.abstol, reltol=opt.reltol,
            callback=cb, tstops=knots, kwargs...)
    return solution_zt
end

# ╔═╡ 33fc9388-a5c9-4b70-91ec-8d4c89bd2026
# normal flow direction
sol1 = solving_migration_time(par1.col.L, par1.col.d, par1.col.df, par1.col.gas, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀, 0.0, par1.sub[1].t₀, par1.opt; t_end=1.0)

# ╔═╡ 415ef581-f292-4209-8087-bbae3bfc4c83
sol1.u[end]

# ╔═╡ 21c81b9a-8542-4228-89c0-f72b9913c504
sol1.retcode # reached t_end = 1.0 s

# ╔═╡ 09dafd05-0549-4945-be77-c79b1402cf67
# reverse flow direction, starting where normal flow stoped after 1s
sol2 = solving_migration_time(par2.col.L, par2.col.d, par2.col.df, par2.col.gas, par2.prog.T_itp, par2.prog.Fpin_itp, par2.prog.pout_itp, par2.sub[1].Tchar, par2.sub[1].θchar, par2.sub[1].ΔCp, par2.sub[1].φ₀, sol1.u[end], par2.sub[1].t₀, par2.opt)

# ╔═╡ 695c5ff8-2a1d-4bf2-9426-2fc04b175148
plot(sol1, label="sol1"); plot!(sol2, label="sol2", legend=:left); xlims!(0.0, 1.0)

# ╔═╡ ed1547b8-fd25-41c3-bbc8-d5df8d428a12
sol2.retcode # reached column inlet x = 0

# ╔═╡ 7712025f-b6f9-479f-8243-ff98c7acebba
sol1_ = solving_migration_time(par1.col.L, par1.col.d, par1.col.df, par1.col.gas, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀, 0.0, par1.sub[1].t₀, par1.opt)

# ╔═╡ d1986e3a-c229-4343-ad38-f8c29c3041db
sol1_.t[end]

# ╔═╡ 8588c9bc-e311-4f07-8d9b-15f5187ee6f2
sol_classic.u[end] == sol1_.t[end]

# ╔═╡ 89708f9b-00da-4abc-922f-afbae8791db3
sol_classic.u[end] - sol1_.t[end] # inside ode-solver tolerance

# ╔═╡ b1ea338f-74be-4faa-a23e-64d866e731c1
sol_classic.u[end] ≈ sol1_.t[end]

# ╔═╡ 63f29ed4-1788-4c9d-b06c-c7d829253712
isapprox(sol_classic.u[end], sol1_.t[end]; atol=opt.abstol, rtol=opt.reltol)

# ╔═╡ 74afc06a-63e8-4cb2-9a1e-81cd20f1dc25
# no flow, pin=pout:
sol3 = solving_migration_time(par3.col.L, par3.col.d, par3.col.df, par3.col.gas, par3.prog.T_itp, par3.prog.Fpin_itp, par3.prog.pout_itp, par3.sub[1].Tchar, par3.sub[1].θchar, par3.sub[1].ΔCp, par3.sub[1].φ₀, sol1.u[end], par3.sub[1].t₀, par3.opt)
# solute stays at position

# ╔═╡ a1a98da5-c246-4c87-864c-f1bd337c38d0
sol3.u[1] - sol3.u[end]

# ╔═╡ 9f8acc89-708b-4fc1-a050-d4f594606702
function solving_bandvariance_time(solution_zt, L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, σ₀, opt; kwargs...)
	z(t) = solution_zt(t)
	p = (L, d, df, gas, T_itp, Fpin_itp, pout_itp, Tchar, θchar, ΔCp, φ₀, Cag, opt)
	f_σ²t(σ²,p,t) = bandode(z(t), t, σ², p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13])
	σ²₀ = σ₀^2
	tspan = (solution_zt.t[begin], solution_zt.t[end])
	knots = program_knots(T_itp, solution_zt.t[begin], solution_zt.t[end])
	
	prob_σ²t = ODEProblem(f_σ²t, σ²₀, tspan, p)
	solution_σ²t = solve(prob_σ²t, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol, tstops=knots; kwargs...)
	return solution_σ²t
end

# ╔═╡ b10ec11a-47dd-486c-86e4-711c58c68dd7
sigma1 = solving_bandvariance_time(sol1, par1.col.L, par1.col.d, par1.col.df, par1.col.gas, par1.prog.T_itp, par1.prog.Fpin_itp, par1.prog.pout_itp, par1.sub[1].Tchar, par1.sub[1].θchar, par1.sub[1].ΔCp, par1.sub[1].φ₀, par1.sub[1].Cag, 0.0, par1.opt)

# ╔═╡ fe96f653-8f5e-4d72-ac02-1f1853002a5b
sqrt(sigma1.u[end])

# ╔═╡ e9f367c3-8a74-4c35-b9b8-29bf0a5c023e
# start with band variance from 1st solution
sigma2 = solving_bandvariance_time(sol2, par2.col.L, par2.col.d, par2.col.df, par2.col.gas, par2.prog.T_itp, par2.prog.Fpin_itp, par2.prog.pout_itp, par2.sub[1].Tchar, par2.sub[1].θchar, par2.sub[1].ΔCp, par2.sub[1].φ₀, par2.sub[1].Cag, sqrt(sigma1.u[end]), par2.opt)

# ╔═╡ 9b9ac87f-3b23-4c66-8f3a-f9323e38b4d5
plot(sigma1, label="sigma1"); plot!(sigma2, label="sigma2")

# ╔═╡ f841aa39-c56d-4cc0-9acd-934845ec5029
# no flow, pin=pout, solute stays at position:
sigma3 = solving_bandvariance_time(sol3, par3.col.L, par3.col.d, par3.col.df, par3.col.gas, par3.prog.T_itp, par3.prog.Fpin_itp, par3.prog.pout_itp, par3.sub[1].Tchar, par3.sub[1].θchar, par3.sub[1].ΔCp, par3.sub[1].φ₀, par3.sub[1].Cag, sqrt(sigma1.u[end]), par3.opt)

# ╔═╡ 47cd2f10-e78f-4323-ae97-26469d82544c
sigma3.u[1] - sigma3.u[end]

# ╔═╡ 7b20d043-0172-4ed5-b663-737699c30c3c
function solving_odesystem_time(
	L, d, df, gas, 
	T_itp, Fpin_itp, pout_itp, 
	Tchar, θchar, ΔCp, φ₀, Cag, 
	t₀, x₀, σ₀, 
	opt; 
	t_end=sum(T_itp.grid[2]), kwargs...
)
	# initial values
	u0 = [x₀; σ₀^2]
	tspan = (t₀, t_end)
	knots = program_knots(T_itp, t₀, t_end)

	p = (
		L, d, df, gas, 
		T_itp, Fpin_itp, pout_itp, 
		Tchar, θchar, ΔCp, φ₀, Cag, 
		opt
	)
    cb = time_domain_exit_callback(L, t₀, p)
	
	prob = ODEProblem(odesystem_time!, u0, tspan, p)
    solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol, callback=cb, tstops=knots; kwargs...)
	return solution
end

# ╔═╡ faa536b2-3b5a-49c4-967c-7215ee7534db
sol_sys1 = solving_odesystem_time(
	par1.col, 
	par1.prog,
	par1.sub[1],
	0.0, 0.0,
	par1.opt;
	t_end=1.0
)

# ╔═╡ cc94b2ae-d76e-4a7f-aa7b-96237e10a2ac
sol_sys1.u[end][1], sol_sys1.u[end][2]

# ╔═╡ 01ed87e8-0b52-4a19-bddc-2bce3c51111a
# reverse flow from 1st solution
sol_sys2 = solving_odesystem_time(
	par2.col, 
	par2.prog,
	par2.sub[1],
	sol_sys1.u[end][1], sqrt(sol_sys1.u[end][2]),
	par2.opt
)

# ╔═╡ dbbed30e-3041-47ba-a9b7-8d835868c331
# no flow after 1st solution, pin=pout, solute stays at position:
sol_sys3 = solving_odesystem_time(
	par3.col, 
	par3.prog,
	par3.sub[1],
	sol_sys1.u[end][1], sqrt(sol_sys1.u[end][2]),
	par3.opt
)

# ╔═╡ 54de9124-42c3-4be1-ac8c-b77ad5ae40dc
sol_sys_1_ = solving_odesystem_time(
	par1.col, 
	par1.prog,
	par1.sub[1],
	0.0, 0.0,
	par1.opt
)

# ╔═╡ 06d90036-ad84-4176-94c9-8843587461df
sol_sys_1_.t[end] # retention time

# ╔═╡ 4924e0e3-0b7f-4cd6-a460-5620afec38df
sol_sys_classic.u[end][1] - sol_sys_1_.t[end]

# ╔═╡ a135841f-7d15-47bb-ba2e-448db14dec3d
sqrt(sol_sys_1_.u[end][2]) # band width in m

# ╔═╡ 86032a6e-f9c1-4337-bb3e-094d4ec13cc3
u_s = uS1(par1.col.L, sol_sys_1_.t[end]) # outlet solute velocity in m/s

# ╔═╡ 13100cd0-e82f-4d8c-862f-615ecb3bc3f9
sqrt(sol_sys_1_.u[end][2]) / u_s # peak width at outlet in s

# ╔═╡ 444af8db-17f8-40ec-9db7-6826e9ee4a82
sqrt(sol_sys_classic.u[end][2]) - sqrt(sol_sys_1_.u[end][2]) / u_s

# ╔═╡ 2c3adb40-fb51-428d-825e-83507242fe78
sol_sys_vac = solving_odesystem_time(
	par_vac.col, 
	par_vac.prog,
	par_vac.sub[1],
	0.0, 0.0,
	par_vac.opt
)

# ╔═╡ 8461f0f5-4a46-41d8-9ec4-6036644b9267
uS_vac(1.0,sol_sys_vac.t[end])

# ╔═╡ 26893e5b-61d2-48b2-b5d4-b12ec700645b
sqrt(sol_sys_vac.u[end][2])/uS_vac(1.0,sol_sys_vac.t[end]) # peak width in s

# ╔═╡ 0ece503d-d7b1-45ac-bc2c-6e0cb9517bc5
sol_sys_new1 = solving_odesystem_time(
	par1.col, 
	par1.prog,
	sub_new,
	par1.opt;
	t_end=1.3
)

# ╔═╡ cacebc8b-9299-4619-a0db-689cf44d00c5
sol_sys_new1.t

# ╔═╡ 8d24472d-0c8c-413f-9788-b428afc0efaf
# reverse flow
begin
	sub_new2 = Substance(sub_new.name, sub_new.CAS, sub_new.Tchar, sub_new.θchar, sub_new.ΔCp, sub_new.φ₀, sub_new.ann, sub_new.Cag, sol_sys_new1.t[end], 0.0, sol_sys_new1.u[end][1], sqrt(sol_sys_new1.u[end][2]))
sol_sys_new2 = solving_odesystem_time(
	par1.col, 
	par2.prog,
	sub_new2,
	par1.opt;
	t_end=100.0
)
end

# ╔═╡ Cell order:
# ╠═423a3ec8-2ed1-46a8-9c2d-57b35d3db73f
# ╠═31e1003e-cc5d-11f1-ac55-13ff6939a21b
# ╠═a1e55ea0-6417-4009-931d-8ce2a2a48ecc
# ╠═09a22481-7dba-460a-b7db-e15bd37607e8
# ╠═56263d3f-dfdc-4863-b16a-decb6ae5c078
# ╠═89923dd0-3ac6-4b73-ae22-4fdf16c5e9b1
# ╠═c6ddea2e-e836-4787-bbe2-6b6892f0366d
# ╠═d508f096-23d8-4275-a7dc-75756b0ebaec
# ╠═8306a4d2-d3ec-44fc-b384-0dfec29502bb
# ╠═cd5b925a-15d9-4470-a9c1-5d6bb27979c3
# ╠═f98a6634-119b-4c57-ba52-336f08c3f354
# ╠═67b656fa-deaa-454a-9954-2c4ca61884d1
# ╠═04ff0803-eb2e-49fb-977d-742093f201f3
# ╠═ebb32727-9603-4fe9-9a06-42ba80a16fff
# ╠═b297bbbf-eb7d-4353-9a5e-e59f30af34d4
# ╠═9cec64ad-3645-4447-bd93-391890bc0422
# ╠═195ed35a-61c9-43ea-807e-a8c7f187f8c0
# ╠═7388cf27-8603-454c-bd9e-92469a7fda02
# ╠═049a20a9-5b86-467a-a752-c022af50ec6c
# ╠═0b10c438-0104-4ba6-9491-9610aa8e5181
# ╠═ddd701f6-b8e9-4cf6-b986-f4f631d456b6
# ╠═7153bf12-e64b-4679-8538-1ec44f151fa4
# ╠═85edb9e5-ee3a-4e36-acb9-2dd036584730
# ╠═74f14b83-3e25-4813-b943-bbf724814f0e
# ╠═8dbae989-f522-456e-a051-6bfe34eebabe
# ╠═38316056-091c-4097-a551-7439f9b4c995
# ╠═d30479bd-8fc2-4990-b1ee-cc5d506e8d7c
# ╠═27e86e63-ed5c-44ce-bad8-74ae2f55c21c
# ╠═de03478f-7af7-40db-8e57-e09552db1667
# ╠═aa5d9359-8771-4c0a-ac98-4dc5e0b34fdc
# ╠═78cae962-2eb3-4bb0-a6a7-62457f16f85d
# ╠═be4cd698-2ba3-4666-bbbd-e91cd9864e47
# ╠═0d5946bd-2fad-4bcd-8bb6-264b4d895a1d
# ╠═0655e3d3-f4cb-4224-8c05-64d7cdb9e3e3
# ╠═07875f26-915e-41b5-8aff-589d891c4f99
# ╠═52ec2ab2-3f37-48d9-be9e-73289ed86ff3
# ╠═9fd00d96-0589-403f-bbf3-2bc6ed7b1e24
# ╠═5d154824-9280-4d7f-9758-58d96fc4136b
# ╠═3262fb55-48b4-407b-a844-df8248afa789
# ╠═7beafde2-7752-4dff-a177-4c2e0cad8a1f
# ╠═768295ef-1658-4ab6-ba85-59945cb8a8dd
# ╠═8a7e4645-9504-415a-a345-cce6ef5c1808
# ╠═f1bd774e-137f-4895-8517-e4f87252649c
# ╠═15f972bb-58d4-49b1-b321-a2f82f322fd4
# ╠═1a8fcf4d-9cf5-41c2-97dc-78c2a66f126a
# ╠═a7bac514-b370-4ed2-8b0e-0cbfd9b3b0bc
# ╠═7b461249-425f-4bf9-ab76-d419ce61ccab
# ╠═a0964a38-0124-47f3-9f12-4f25a25007f6
# ╠═dbcfa37e-02f8-4b1f-9871-6936e53052be
# ╠═b8282421-962d-4ee2-b65e-8b6f260bde1e
# ╠═ff76d100-9c5b-4068-93e3-62790f45549a
# ╠═e839ccda-04cc-4f2b-9a60-248a230b17d2
# ╠═62533e44-c765-4a54-932d-0a6d6b2024d6
# ╠═d6c17996-1ead-4281-bc7e-7629b0afb745
# ╠═fb02ba01-0e4f-4e1a-8bc1-5e805a0e1f40
# ╠═81a48945-ea78-4ac9-b128-c8426cd7c071
# ╠═5ad305c4-dfd5-40fd-b422-cda482ab7acf
# ╠═dcbb2198-66eb-48d1-a4a1-54375de2bf98
# ╠═0654d56f-2458-450f-995d-257c86abe2a3
# ╠═c08c8663-e257-4228-b358-c5f40320efb6
# ╠═fb9f99af-183b-45db-8868-5b2a199f8916
# ╠═4844421b-c7a7-4029-8dd5-dba9ab829b88
# ╠═570e1699-86f9-43d9-812a-f78bf1402b54
# ╠═a43dc802-b93a-4a17-9135-67cb5a800297
# ╠═194dc600-bb61-4040-8ae0-78c2f92ff09b
# ╠═3cd03859-8aaf-4ad1-a1f1-84534d42a385
# ╠═b7fbab24-fb3d-42fd-9c53-9e1d510ba69c
# ╠═2a4c8410-b5ec-4a80-95a6-22e0d264aa99
# ╠═04199e01-eac5-4b01-8a4a-504fe6b415ad
# ╠═144dc9b7-d0ab-4f80-846c-165040f299f3
# ╠═0236dc1c-e1ad-480c-8e69-c087bd1085cf
# ╠═85f94ab5-2135-4641-ab46-b007214b91ee
# ╠═d78e12a1-b46d-4aed-80fb-9c0c4deb3983
# ╠═e51628de-5a58-4ec2-8d37-adf435b276dd
# ╠═6e78cdef-7f2f-4ed4-abce-0964905ade80
# ╠═04f05dfc-4ea2-4ffe-a22c-6f06aa81e367
# ╠═302691e9-01d9-44cc-9da6-a757fbe13012
# ╠═377cc5c2-280f-4d0e-b8d5-f467a34edc99
# ╠═00d98412-72e6-4704-9ca0-c22c1a7937bb
# ╠═404ef427-7785-4649-bc1b-143c432c27de
# ╠═28421f8d-a107-42ac-b222-da21d5fec223
# ╠═b48d5a49-ae01-41da-8983-ba522ea7b811
# ╠═52f0f5c6-a3f2-4c71-8f8b-8fc92f7f9694
# ╠═0cf79001-43d7-432f-ba24-cbfb6ba3b93e
# ╠═9480b569-3a5b-4f04-9acf-b91235120fc7
# ╠═d27b6674-6857-441f-ab7a-0ea13ae83690
# ╠═5f134571-2c33-4a2c-99e7-06050544c6ee
# ╠═9b772052-82ad-4fcc-b35e-7c9da66ce6cc
# ╠═8ca13a43-0b0a-43a5-a9d6-4e61a54299ab
# ╠═a6d98683-98e5-4074-8ca6-c8eb89fd17f3
# ╠═d24a67e8-e01b-4b01-aea8-2aed0f5cbf13
# ╠═2cf1c3f4-b093-4cbb-a2a4-55d429ec841e
# ╠═a95a9a9f-e1a2-4c2a-aebc-ff100ca10c9e
# ╠═25a4a2c3-ed68-4d4a-9432-298d3a5fe27e
# ╠═a8471106-8e88-4abc-8841-2ab791a72f77
# ╠═1170bcc9-ac3d-44af-9c67-cbd9b85a7e3e
# ╠═a3e1934a-8b28-404b-9a4f-c910bddce852
# ╠═3fc58fe7-8387-4952-ab50-d49640bc1222
# ╠═f18f1d88-dbbb-4528-a803-607d60adf84b
# ╠═a4d85b2e-197f-41f4-a602-e0bb83204cdf
# ╠═c08f9403-61bb-4478-81ec-8cd6d73cf16a
# ╠═33fc9388-a5c9-4b70-91ec-8d4c89bd2026
# ╠═415ef581-f292-4209-8087-bbae3bfc4c83
# ╠═695c5ff8-2a1d-4bf2-9426-2fc04b175148
# ╠═21c81b9a-8542-4228-89c0-f72b9913c504
# ╠═09dafd05-0549-4945-be77-c79b1402cf67
# ╠═ed1547b8-fd25-41c3-bbc8-d5df8d428a12
# ╠═00bd501f-ca26-47ae-a39c-b48b12121b99
# ╠═7712025f-b6f9-479f-8243-ff98c7acebba
# ╠═5e7746c3-c15e-427a-a73a-fbc0120e51a2
# ╠═d1986e3a-c229-4343-ad38-f8c29c3041db
# ╠═8588c9bc-e311-4f07-8d9b-15f5187ee6f2
# ╠═89708f9b-00da-4abc-922f-afbae8791db3
# ╠═b1ea338f-74be-4faa-a23e-64d866e731c1
# ╠═63f29ed4-1788-4c9d-b06c-c7d829253712
# ╠═74afc06a-63e8-4cb2-9a1e-81cd20f1dc25
# ╠═a1a98da5-c246-4c87-864c-f1bd337c38d0
# ╠═a1f284b5-5d46-44aa-b954-c768040076b6
# ╠═3579dae0-c0ed-4775-a134-acde3090c50e
# ╠═9f8acc89-708b-4fc1-a050-d4f594606702
# ╠═71fd80e1-6784-433a-949e-a507945efdc6
# ╠═1e0aa0b0-9f68-4458-9644-1403aa036a15
# ╠═abd96ffd-faf4-4d34-9790-42014aefddd8
# ╠═9f41ed53-5a28-4862-b172-fdd56da6612b
# ╠═b10ec11a-47dd-486c-86e4-711c58c68dd7
# ╠═9b9ac87f-3b23-4c66-8f3a-f9323e38b4d5
# ╠═fe96f653-8f5e-4d72-ac02-1f1853002a5b
# ╠═e9f367c3-8a74-4c35-b9b8-29bf0a5c023e
# ╠═f841aa39-c56d-4cc0-9acd-934845ec5029
# ╠═47cd2f10-e78f-4323-ae97-26469d82544c
# ╠═54eeb56a-b6f8-4bcd-a477-46cd61a0ec0b
# ╠═7b20d043-0172-4ed5-b663-737699c30c3c
# ╠═f5db1295-221b-4cf2-8d8d-f076610c22c3
# ╠═5e7a2e30-dc78-4a8e-9284-e74e3414e898
# ╠═efbf0fd3-8fb5-4292-851e-14d68472ef75
# ╠═c14f63c8-3708-4d83-bbcf-596b2d88057d
# ╠═faa536b2-3b5a-49c4-967c-7215ee7534db
# ╠═cc94b2ae-d76e-4a7f-aa7b-96237e10a2ac
# ╠═01ed87e8-0b52-4a19-bddc-2bce3c51111a
# ╠═dbbed30e-3041-47ba-a9b7-8d835868c331
# ╠═8617ab5b-8657-484e-9e10-41055db08033
# ╠═54de9124-42c3-4be1-ac8c-b77ad5ae40dc
# ╠═7be3558a-deed-4f05-82a8-b50906b73712
# ╠═06d90036-ad84-4176-94c9-8843587461df
# ╠═4924e0e3-0b7f-4cd6-a460-5620afec38df
# ╠═76c707f2-7b44-4058-8a3a-3af98dee3bb2
# ╠═a135841f-7d15-47bb-ba2e-448db14dec3d
# ╠═86032a6e-f9c1-4337-bb3e-094d4ec13cc3
# ╠═13100cd0-e82f-4d8c-862f-615ecb3bc3f9
# ╠═444af8db-17f8-40ec-9db7-6826e9ee4a82
# ╠═ec119d1b-0bca-4ae6-aa97-5839ae5f07dd
# ╠═95ccf180-c5ed-4ae8-a4aa-71e7f79a7f60
# ╠═2c3adb40-fb51-428d-825e-83507242fe78
# ╠═49abc3ae-ff8a-4869-8f9d-50c3559c40cc
# ╠═0a3ed9d8-40d8-47c0-87ac-84e4be66be21
# ╠═8461f0f5-4a46-41d8-9ec4-6036644b9267
# ╠═26893e5b-61d2-48b2-b5d4-b12ec700645b
# ╠═c214a16a-4055-4bca-91ca-9098477ea847
# ╠═494b6f1a-13a2-420a-958b-a120237e6092
# ╠═7cdaa4a9-eed0-4fc6-959d-755e0602fa60
# ╠═64b9dd0c-7721-4df1-8a92-6b3a8541887b
# ╠═e6e023a2-9d18-411b-8fb2-a3d682d56675
# ╠═3832d314-4e48-4be7-b00b-9ad588c12c54
# ╠═0ece503d-d7b1-45ac-bc2c-6e0cb9517bc5
# ╠═cacebc8b-9299-4619-a0db-689cf44d00c5
# ╠═b2905bc6-00f7-4bd6-8e44-4896ea14142f
# ╠═8d24472d-0c8c-413f-9788-b428afc0efaf
# ╠═100ef7fc-3162-4872-af25-fa132ae94ebc
# ╠═562a83c8-2373-4bce-ac89-189174573d95
# ╠═a75817af-3f1a-4e55-bd1c-f84a3e8a26f6
