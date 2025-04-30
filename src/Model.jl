#---Begin-Functions-of-the-physical-model---
"""
    pressure(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")

Calculate the pressure at position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `Fpin_itp`: Interpolated (linear) flow `F(t)` resp. inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m. Can be a function of position `x`.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
    or with a gradient (`ng = false`, eq. 1).
* `vis`: used model for viscosity "Blumberg" or "HP"
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``p(x,t) =
\\sqrt(p_{in}(t)^2-\\frac{κ(x,t)}{κ_L(t)}\\left(p_{in}^2-p_{out}^2\\right))``
Eq. 1

``p(x,t) =
\\sqrt(p_{in}(t)^2-\\frac{x}{L}\\left(p_{in}^2-p_{out}^2\\right))`` Eq. 2

with ``κ(x,t)`` the flow restriction up to position `x` at time `t` and
``κ_L(t) = κ(x=L,t)`` the flow restriction of the whole column at
time `t`.

See also: [`flow_restriction`](@ref)
"""
function pressure(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")
    if ng==true
        if control == "Pressure"
            pin_itp = Fpin_itp
            pp = real(sqrt(Complex(pin_itp(t)^2 - x/L*(pin_itp(t)^2-pout_itp(t)^2))))
        elseif control == "Flow"
            F_itp = Fpin_itp
            if isa(d, Function) 
                pp = real(sqrt(Complex(pout_itp(t)^2 + 256/π * pn/Tn * viscosity(x, t, T_itp, gas, vis=vis)* T_itp(x,t)/d(x)^4 * F_itp(t) * (L - x))))
            else
                pp = real(sqrt(Complex(pout_itp(t)^2 + 256/π * pn/Tn * viscosity(x, t, T_itp, gas, vis=vis)* T_itp(x,t)/d^4 * F_itp(t) * (L - x))))
            end
        end
    else
        if control == "Pressure"
            pin_itp = Fpin_itp
            pp = real(sqrt(Complex(pin_itp(t)^2 - flow_restriction(x, t, T_itp, d, gas, vis=vis)/flow_restriction(L, t, T_itp, d, gas, vis=vis)*(pin_itp(t)^2-pout_itp(t)^2))))
        elseif control == "Flow"
            F_itp = Fpin_itp
            κL = flow_restriction(L, t, T_itp, d, gas, vis=vis)
            κx = flow_restriction(x, t, T_itp, d, gas, vis=vis)
            pp = real(sqrt(Complex(pout_itp(t)^2 + 256/π * pn/Tn * F_itp(t) * (κL - κx))))
        end
    end
    return pp
end

function pressure(x, t, par::Parameters)
    pp = pressure(x, t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
    return pp
end

"""
    flow_restriction(x, t, T_itp, d, gas; ng=false, vis="Blumberg")

Calculate the flow restriction ``κ`` up to position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `d`: Diameter of the GC column, in m. Can be a function of position `x`.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
    or with a gradient (`ng = false`, eq. 1).
* `vis`: used model for viscosity "Blumberg" or "HP"

``κ(x,t) = \\int_0^x \\frac{η(y,t) T(y,t)}{d(y)^4}dy``
Eq. 1

``κ(x,t) = \\frac{η(t) T(t) x}{d^4}`` Eq. 2

with ``η(x,t)`` the viscosity of the mobile phase gas.

See also: [`viscosity`](@ref)
"""
function flow_restriction(x, t, T_itp, d, gas; ng=false, vis="Blumberg")
    if ng==true
        κ = x*GasChromatographySimulator.viscosity(x, t, T_itp, gas, vis=vis)*T_itp(x, t)*d(x)^-4
    else
		p = (t, )# d) if `d` is function of `x` it cannot be part of `p` 
        f(y, p) = GasChromatographySimulator.viscosity(y, p[1], T_itp, gas, vis=vis)*T_itp(y, p[1])*d(y)^-4
        domain = (zero(x), x)  # Use zero(x) to get a zero of the same type as x
        prob = IntegralProblem(f, domain, p)
        κ = solve(prob, GasChromatographySimulator.QuadGKJL(), reltol=1e-3, abstol=1e-3)[1]
    end
    return κ
end

function flow_restriction(x, t, T_itp, d::Number, gas; ng=false, vis="Blumberg")
    if ng==true
        κ = x*GasChromatographySimulator.viscosity(x, t, T_itp, gas, vis=vis)*T_itp(x, t)*d^-4
    else
		p = (t, d)
        f(y, p) = GasChromatographySimulator.viscosity(y, p[1], T_itp, gas, vis=vis)*T_itp(y, p[1])*p[2]^-4
        domain = (zero(x), x)  # Use zero(x) to get a zero of the same type as x
        prob = IntegralProblem(f, domain, p)
        κ = solve(prob, GasChromatographySimulator.QuadGKJL(), reltol=1e-3, abstol=1e-3)[1]
    end
    return κ
end

function flow_restriction(x, t, par::Parameters)
    κ = flow_restriction(x, t, par.prog.T_itp, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis)
    return κ
end

"""
    viscosity(x, t, T_itp, gas; vis="Blumberg")

Calculate the (dynamic) viscosity of the mobile phase gas at position `x`
at time `t` in Pa s.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `gas`: Name of the mobile phase gas.
* `vis`: used model

`vis = "HP"`

Simple model used in the HP Flow calculator

``η(x,t) = C_1 * \\left(T(x,t) + T_{st}\\right) + C_2`

`vis = "Blumberg"`

``η(x,t) = η_{st}\\left(\\frac{T(x,t)}{T_{st}}\\right)^{(ξ_0 + ξ_1 \\frac{T(x,t)-T_{st}}{T_{st}})}`` 

with ``η_{st}``, ``ξ_0`` and ``ξ_1`` parameters dependent on the
mobile phase gas [1].

[1] Blumberg, Leonid M., Temperature-Programmed Gas Chromatography,
Wiley-VCH, 2010.
"""
function viscosity(x, t, T_itp, gas; vis="Blumberg")
    if vis == "Blumberg"
        if gas=="He"
            ηst = 18.63e-6
            ξ₀ = 0.6958
            ξ₁ = -0.0071
        elseif gas=="H2"
            ηst = 8.382e-6
            ξ₀ = 0.6892
            ξ₁ = 0.005
        elseif gas=="N2"
            ηst = 16.62e-6
            ξ₀ = 0.7665
            ξ₁ = -0.0378
        elseif gas=="Ar"
            ηst = 21.04e-6
            ξ₀ = 0.8131
            ξ₁ = -0.0426
        else
            error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
        end
        T = T_itp(x, t)
        η = ηst*(T/Tst)^(ξ₀ + ξ₁*(T-Tst)/Tst)
    elseif vis == "HP"
        if gas=="He"
            C₁ = 4.28e-8
            C₂ = 6.968e-6
        elseif gas=="H2"
            C₁ = 3.5e-8
            C₂ = 7.994e-6
        elseif gas=="N2"
            C₁ = 1.83e-8
            C₂ = 4.416e-6
        else
            error("Unknown selection of gas. Choose one of these: He, H2 or N2.")
        end
        T = T_itp(x, t)
        η = C₁*T + C₂
    else
        error("Unknown selection for the viscosity model. Choose one of these options: 'Blumberg' or 'HP'.")
    end
    return η
end

function viscosity(x, t, par::Parameters)
    η = viscosity(x, t, par.prog.T_itp, par.col.gas; vis=par.opt.vis)
    return η
end

"""
    viscosity(T, gas; vis="Blumberg")

Calculate the (dynamic) viscosity of the mobile phase gas at temperature `T` in Pa s.

# Arguments
* `T`: Temperature in K.
* `gas`: Name of the mobile phase gas.
* `vis`: used model for viscosity "Blumberg" or "HP".

`vis = "HP"`

Simple model used in the HP Flow calculator

``η(x,t) = C_1 * \\left(T(x,t) + T_{st}\\right) + C_2`

`vis = "Blumberg"``

``η(x,t) = η_{st}\\left(\\frac{T)}{T_{st}}\right)^{(ξ_0 + ξ_1 \\frac{T-T_{st}}{T_{st}})}`` 

with ``η_{st}``, ``ξ_0`` and ``ξ_1`` parameters dependent on the
mobile phase gas [1].

[1] Blumberg, Leonid M., Temperature-Programmed Gas Chromatography,
Wiley-VCH, 2010.
"""
function viscosity(T, gas::String; vis="Blumberg")
    if vis == "Blumberg"
        if gas=="He"
            ηst = 18.63e-6
            ξ₀ = 0.6958
            ξ₁ = -0.0071
        elseif gas=="H2"
            ηst = 8.382e-6
            ξ₀ = 0.6892
            ξ₁ = 0.005
        elseif gas=="N2"
            ηst = 16.62e-6
            ξ₀ = 0.7665
            ξ₁ = -0.0378
        elseif gas=="Ar"
            ηst = 21.04e-6
            ξ₀ = 0.8131
            ξ₁ = -0.0426
        else
            error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
        end
        η = ηst*(T/Tst)^(ξ₀ + ξ₁*(T-Tst)/Tst)
    elseif vis == "HP"
        if gas=="He"
            C₁ = 4.28e-8
            C₂ = 6.968e-6
        elseif gas=="H2"
            C₁ = 3.5e-8
            C₂ = 7.994e-6
        elseif gas=="N2"
            C₁ = 1.83e-8
            C₂ = 4.416e-6
        else
            error("Unknown selection of gas. Choose one of these: He, H2 or N2.")
        end
        η = C₁*T + C₂
    else
        error("Unknown selection for the viscosity model. Choose one of these options: 'Blumberg' or 'HP'.")
    end
    return η
end

"""
    inlet_pressure(t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")

Calculate the inlet pressure for a given column, temperature, flow and outlet pressure
"""
function inlet_pressure(t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")
    if control == "Pressure"
        pin = Fpin_itp(t)
    elseif control == "Flow"
        F_itp = Fpin_itp
        if ng == true 
            pin = real(sqrt(Complex(pout_itp(t)^2 + 256/π * pn/Tn * viscosity(0.0, t, T_itp, gas; vis="Blumberg") * T_itp(0.0, t) * L / d(0.0)^4 * F_itp(t))))
        else
            κL = flow_restriction(L, t, T_itp, d, gas, vis=vis)
            pin = real(sqrt(Complex(pout_itp(t)^2 + 256/π * pn/Tn * κL * F_itp(t))))
        end
    end
    return pin
end

function inlet_pressure(t, T_itp, Fpin_itp, pout_itp, L, d::Number, gas; ng=false, vis="Blumberg", control="Pressure")
    if control == "Pressure"
        pin = Fpin_itp(t)
    elseif control == "Flow"
        F_itp = Fpin_itp
        if ng == true 
            pin = real(sqrt(Complex(pout_itp(t)^2 + 256/π * pn/Tn * viscosity(0.0, t, T_itp, gas; vis="Blumberg") * T_itp(0.0, t) * L / d^4 * F_itp(t))))
        else
            κL = flow_restriction(L, t, T_itp, d, gas, vis=vis)
            pin = real(sqrt(Complex(pout_itp(t)^2 + 256/π * pn/Tn * κL * F_itp(t))))
        end
    end
    return pin
end

function inlet_pressure(t, par::Parameters)
    pin = inlet_pressure(t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
    return pin
end

"""
    holdup_time(T, Fpin, pout, L, d, gas; vis="Blumberg", control="Pressure")

Calculate the hold-up time in s without a gradient.

# Arguments
* `T`: Temperature in K.
* `Fpin`: Flow in m³/s resp. inlet pressure in Pa(a).
* `pout`: Outlet pressure in Pa(g).
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``t_M = \\frac{128}{3}\\frac{L^2}{d^2}η\\frac{p_{in}^3-p_{out}^3}{(p_{in}^2-p_{out}^2)^2}``
"""
function holdup_time(T::Float64, Fpin::Float64, pout::Float64, L::Float64, d::Float64, gas::String; vis="Blumberg", control="Pressure")
    # hold-up time at temperature T (non-gradient)
	η = viscosity(T, gas; vis=vis)
    if control == "Pressure"
        pin = Fpin
    elseif control == "Flow"
        F = Fpin
        pin = real(sqrt(Complex(pout^2 + 256/π * pn/Tn * viscosity(T, gas; vis="Blumberg") * T * L / d^4 * F)))
    end
	tM = 128/3*L^2/d^2*η*(pin^3-pout^3)/(pin^2-pout^2)^2
	return tM
end

"""
    holdup_time(t, T_itp, pin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")

Calculate the hold-up time in s at time `t` with a gradient.

# Arguments
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `pin_itp`: Interpolated (linear) inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
    or with a gradient (`ng = false`, eq. 1).
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``t_M(t) = 64\\frac{κ_L(t)}{p_{in}(t)^2-p_{out}(t)^2} \\int_0^L
d(y)^2\\frac{p(y,t)}{T(y,t)}dy`` Eq. 1

``t_M(t) =
\\frac{128}{3}\\frac{L^2}{d^2}η(t)\\frac{p_{in}(t)^3-p_{out}(t)^3}{(p_{in}(t)^2-p_{out}(t)^2)^2}``
Eq. 2
"""
function holdup_time(t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")
    # hold-up time at time t in a temperature program with potential thermal gradient
    if control == "Pressure"
        #pin(t) = Fpin_itp(t)
        if ng==true
            η = GasChromatographySimulator.viscosity(L, t, T_itp, gas; vis=vis)
            tM = 128/3*L^2/d(L)^2*η*(Fpin_itp(t)^3-pout_itp(t)^3)/(Fpin_itp(t)^2-pout_itp(t)^2)^2
        else
            κL = flow_restriction(L, t, T_itp, d, gas; ng=false, vis=vis)
            f_p(y, p) = d(y)^2*pressure(y, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis=vis, control="Pressure")/T_itp(y, t)
            domain = (0.0, L)
            prob_p = IntegralProblem(f_p, domain)
            integral = solve(prob_p, QuadGKJL(), reltol=1e-3, abstol=1e-3)[1]
            tM = 64*κL/(Fpin_itp(t)^2-pout_itp(t)^2)*integral
        end
    elseif control == "Flow"
        if ng==true
            pin = GasChromatographySimulator.inlet_pressure(t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=true, vis=vis, control="Flow")
            η = GasChromatographySimulator.viscosity(L, t, T_itp, gas; vis=vis)
            tM = 128/3*L^2/d(L)^2*η*(pin^3-pout_itp(t)^3)/(pin^2-pout_itp(t)^2)^2
        else
            f_F(y, p) = d(y)^2*pressure(y, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis=vis, control="Flow")/T_itp(y, t)
            domain = (0.0, L)
            prob_F = IntegralProblem(f_F, domain)
            integral = solve(prob_F, QuadGKJL(), reltol=1e-3, abstol=1e-3)[1]
            tM = π/4 * Tn/pn * integral/Fpin_itp(t)
        end
    end
    
    return tM
end

function holdup_time(t, T_itp, Fpin_itp, pout_itp, L, d::Number, gas; ng=false, vis="Blumberg", control="Pressure")
    # hold-up time at time t in a temperature program with potential thermal gradient
    if control == "Pressure"
        #pin(t) = Fpin_itp(t)
        if ng==true
            η = GasChromatographySimulator.viscosity(L, t, T_itp, gas; vis=vis)
            tM = 128/3*L^2/d^2*η*(Fpin_itp(t)^3-pout_itp(t)^3)/(Fpin_itp(t)^2-pout_itp(t)^2)^2
        else
            κL = flow_restriction(L, t, T_itp, d, gas; ng=false, vis=vis)
            f_p(y, p) = d^2*pressure(y, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis=vis, control="Pressure")/T_itp(y, t)
            domain = (0.0, L)
            prob_p = IntegralProblem(f_p, domain)
            integral = solve(prob_p, QuadGKJL(), reltol=1e-3, abstol=1e-3)[1]
            tM = 64*κL/(Fpin_itp(t)^2-pout_itp(t)^2)*integral
        end
    elseif control == "Flow"
        if ng==true
            pin = GasChromatographySimulator.inlet_pressure(t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=true, vis=vis, control="Flow")
            η = GasChromatographySimulator.viscosity(L, t, T_itp, gas; vis=vis)
            tM = 128/3*L^2/d^2*η*(pin^3-pout_itp(t)^3)/(pin^2-pout_itp(t)^2)^2
        else
            f_F(y, p) = d^2*pressure(y, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis=vis, control="Flow")/T_itp(y, t)
            domain = (0.0, L)
            prob_F = IntegralProblem(f_F, domain)
            integral = solve(prob_F, QuadGKJL(), reltol=1e-3, abstol=1e-3)[1]
            tM = π/4 * Tn/pn * integral/Fpin_itp(t)
        end
    end
    
    return tM
end

function holdup_time(t, par::Parameters)
    tM = holdup_time(t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
    return tM
end

"""
    flow(T, Fpin, pout, L, d, gas; vis="Blumberg", control="Pressure")

Calculate the normalized flow through the GC column in m³/s without a gradient.

# Arguments
* `T`: Temperature in K.
* `Fpin`: Flow in m³/s resp. inlet pressure in Pa(a).
* `pout`: Outlet pressure in Pa(g).
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``F =
\\frac{π}{256}\\frac{T_n}{p_n}\\frac{d^4}{L}\\frac{p_{in}^2-p_{out}^2}{η T}``

with ``T_n`` the normalized temperature (``T_n=(25 + 273.15)``K), ``p_n``
the normalized pressure (``p_n = 101300`` Pa(a)) and ``η`` the viscosity
the mobile phase gas at temperature ``T``.
"""
function flow(T::Float64, Fpin::Float64, pout::Float64, L::Float64, d::Float64, gas::String; vis="Blumberg", control="Pressure")
	# normalized Flow at temperature T (non-gradient)
    if control == "Pressure"
        pin = Fpin
	    η = viscosity(T, gas; vis=vis)
	    F = π/256 * Tn/pn * d^4/L * (pin^2-pout^2)/(η*T)
    elseif control == "Flow"
        F = Fpin
    end
	return F
end

"""
    flow(t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")

Calculate the normalized flow through the GC column in m³/s at time `t`.

# Arguments
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `Fpin_itp`: Interpolated (linear) Flow F(t) resp. inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`,
    eq. 2)
    or with a gradient (`ng = false`, eq. 1).
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``F(t) =
\\frac{π}{256}\\frac{T_n}{p_n}\\frac{p_{in}(t)^2-p_{out}(t)^2}{κ_L(t)}``
Eq. 1

``F(t) =
\\frac{π}{256}\\frac{T_n}{p_n}\\frac{d^4}{L}\\frac{p_{in}(t)^2-p_{out}(t)^2}{η(t)
T(t)}``
Eq. 2

with ``T_n`` the normalized temperature (``T_n=(25 + 273.15)``K), ``p_n``
the normalized pressure (``p_n = 101300`` Pa(a)), ``κ_L`` the flow
restriction of the column and ``η`` the viscosity
the mobile phase gas at temperature ``T``.
"""
function flow(t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")
	# normalized Flow at time t in a temperature program with potential thermal
	# gradient
    # TODO: test for gradient in d(x)
    if control == "Pressure"
        pin_itp = Fpin_itp
        if ng==true
            η = GasChromatographySimulator.viscosity(L, t, T_itp, gas, vis=vis)
            F = π/256 * Tn/pn * d(L)^4/L * (pin_itp(t)^2-pout_itp(t)^2)/(η*T_itp(L,t))
        else
            κL = flow_restriction(L, t, T_itp, d, gas; ng=false, vis=vis)
            F = π/256 * Tn/pn * (pin_itp(t)^2-pout_itp(t)^2)/κL
        end
    elseif control == "Flow"
        F = Fpin_itp(t)
    end
	return F
end

function flow(t, T_itp, Fpin_itp, pout_itp, L, d::Number, gas; ng=false, vis="Blumberg", control="Pressure")
	# normalized Flow at time t in a temperature program with potential thermal
	# gradient
    # TODO: test for gradient in d(x)
    if control == "Pressure"
        pin_itp = Fpin_itp
        if ng==true
            η = GasChromatographySimulator.viscosity(L, t, T_itp, gas, vis=vis)
            F = π/256 * Tn/pn * d^4/L * (pin_itp(t)^2-pout_itp(t)^2)/(η*T_itp(L,t))
        else
            κL = flow_restriction(L, t, T_itp, d, gas; ng=false, vis=vis)
            F = π/256 * Tn/pn * (pin_itp(t)^2-pout_itp(t)^2)/κL
        end
    elseif control == "Flow"
        F = Fpin_itp(t)
    end
	return F
end

function flow(t, par::Parameters)
    F = flow(t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
    return F
end

"""
    mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")

Calculate the residency (the inverse velocity) of the mobile phase at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `Fpin_itp`: Interpolated (linear) Flow F(t) resp. inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``r_M(x,t) = 64 \\frac{d^2 κ_L}{T(x,t)}\\frac{p(x,t)}{p_{in}^2-p_{out}^2}``

with ``T_n`` the normalized temperature (``T_n=(25 + 273.15)``K), ``p_n``
the normalized pressure (``p_n = 101300`` Pa(a)), ``κ_L`` the flow
restriction of the column and ``p(x,t)`` the local pressure.

See also: [`pressure`](@ref), [`flow_restriction`](@ref)
"""
function mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=false, vis="Blumberg", control="Pressure")
    if control == "Pressure"
        pin_itp = Fpin_itp
        pp = pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control="Pressure")
        κL = flow_restriction(L, t, T_itp, d, gas; ng=ng, vis=vis)
        rM = 64*(pp*(d(x))^2)/T_itp(x, t)*κL/(pin_itp(t)^2-pout_itp(t)^2)
    elseif control == "Flow"
        F_itp = Fpin_itp
        pp = pressure(x, t, T_itp, F_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control="Flow")
        rM = π/4 * Tn/pn * d(x)^2/F_itp(t) * pp/T_itp(x,t)
    end
    return rM
end

function mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d::Number, gas; ng=false, vis="Blumberg", control="Pressure")
    if control == "Pressure"
        pin_itp = Fpin_itp
        pp = pressure(x, t, T_itp, pin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control="Pressure")
        κL = flow_restriction(L, t, T_itp, d, gas; ng=ng, vis=vis)
        rM = 64*(pp*d^2)/T_itp(x, t)*κL/(pin_itp(t)^2-pout_itp(t)^2)
    elseif control == "Flow"
        F_itp = Fpin_itp
        pp = pressure(x, t, T_itp, F_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control="Flow")
        rM = π/4 * Tn/pn * d^2/F_itp(t) * pp/T_itp(x,t)
    end
    return rM
end

function mobile_phase_residency(x, t, par::Parameters)
    rM = mobile_phase_residency(x, t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
    return rM
end

"""
    residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp,  φ₀; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)

Calculate the residency (the inverse velocity) of the solute at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `Fpin_itp`: Interpolated (linear) Flow F(t) resp. inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `df`: Film thickness of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `Tchar`: Characteristic temperature of the solute, in K.
* `θchar`: Characteristic parameters of the solute, in °C.
* `ΔCp`: Change of the isobaric heat capacity of the solute moving from the mobile to the
stationary phase, in J mol⁻¹ K⁻¹.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the
thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``r(x,t) = r_M(x,t) \\left(1+k(x,t)\\right)``

with ``r_M`` the residency of the mobile phase and ``k(x,t)`` the retention
factor of the solute on the stationary phase.

See also: [`mobile_phase_residency`](@ref), [`retention_factor`](@ref)
"""
function residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    r = mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)*(1 + retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀; k_th=k_th))
    return r
end

function residency(x, t, col::Column, prog::Program, sub::Substance, opt::Options)
    r = residency(x, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
    return r
end

"""
    retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀)

Calculate the retention factor of the solute in the stationary phase at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `d`: Diameter of the GC column, in m.
* `df`: Film thickness of the GC column, in m.
* `Tchar`: Characteristic temperature of the solute, in K.
* `θchar`: Characteristic parameters of the solute, in °C.
* `ΔCp`: Change of the isobaric heat capacity of the solute moving from the mobile to the
stationary phase, in J mol⁻¹ K⁻¹.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the
thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.

``k(x,t) = \\frac{φ}{φ₀}
\\exp{\\left((\\frac{ΔC_p}{R}+\\frac{T_{char}}{θ_{char}})(\\frac{T_{char}}{T}+-1)
    \\frac{ΔC_p}{R}\\ln{(\\frac{T}{T_{char}})}\\right)}``

with ``R`` the molar gas constant and ``φ`` the dimensionless film thickness
of the simulated GC Column (``φ = d_f/d``).

**TODO**: add option for the retention model ('ABC', 'K-centric')
"""
function retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀; k_th=1e12)
    # this version of the function, where every parameter is
    # given to the function separatly seems to be the fastest
    # version
    # for now only the ideal thermodynamic model
    if Tchar == 0.0 && θchar == 0.0 && ΔCp == 0.0 # non-retained solute
        k = 0.0
    else
        T = T_itp(x, t)
        φ = df(x)/d(x)
        C = ΔCp/R
        lnk₀ = (C + Tchar/θchar) * (Tchar/T - 1) + C*log(T/Tchar)
        k = φ/φ₀*exp(lnk₀)
    end
    if k > k_th
        k = k_th
    end
    return k
end

function retention_factor(x, t, T_itp, d::Number, df::Number, Tchar, θchar, ΔCp, φ₀; k_th=1e12)
    # this version of the function, where every parameter is
    # given to the function separatly seems to be the fastest
    # version
    # for now only the ideal thermodynamic model
    if Tchar == 0.0 && θchar == 0.0 && ΔCp == 0.0 # non-retained solute
        k = 0.0
    else
        T = T_itp(x, t)
        φ = df/d
        C = ΔCp/R
        lnk₀ = (C + Tchar/θchar) * (Tchar/T - 1) + C*log(T/Tchar)
        k = φ/φ₀*exp(lnk₀)
    end
    if k > k_th
        k = k_th
    end
    return k
end

function retention_factor(x, t, T_itp, d, df::Number, Tchar, θchar, ΔCp, φ₀; k_th=1e12)
    # this version of the function, where every parameter is
    # given to the function separatly seems to be the fastest
    # version
    # for now only the ideal thermodynamic model
    if Tchar == 0.0 && θchar == 0.0 && ΔCp == 0.0 # non-retained solute
        k = 0.0
    else
        T = T_itp(x, t)
        φ = df/d(x)
        C = ΔCp/R
        lnk₀ = (C + Tchar/θchar) * (Tchar/T - 1) + C*log(T/Tchar)
        k = φ/φ₀*exp(lnk₀)
    end
    if k > k_th
        k = k_th
    end
    return k
end

function retention_factor(x, t, T_itp, d::Number, df, Tchar, θchar, ΔCp, φ₀; k_th=1e12)
    # this version of the function, where every parameter is
    # given to the function separatly seems to be the fastest
    # version
    # for now only the ideal thermodynamic model
    if Tchar == 0.0 && θchar == 0.0 && ΔCp == 0.0 # non-retained solute
        k = 0.0
    else
        T = T_itp(x, t)
        φ = df(x)/d
        C = ΔCp/R
        lnk₀ = (C + Tchar/θchar) * (Tchar/T - 1) + C*log(T/Tchar)
        k = φ/φ₀*exp(lnk₀)
    end
    if k > k_th
        k = k_th
    end
    return k
end

function retention_factor(x, t, col::Column, prog::Program, sub::Substance, opt::Options)
    k = retention_factor(x, t, prog.T_itp, col.d, col.df, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀; k_th=opt.k_th)
    return k
end

"""
    plate_height(x, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)

Calculate the plate height of the solute at position `x` at time `t`
according to the Golay equation.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `Fpin_itp`: Interpolated (linear) Flow F(t) resp. inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `df`: Film thickness of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `Tchar`: Characteristic temperature of the solute, in K.
* `θchar`: Characteristic parameters of the solute, in °C.
* `ΔCp`: Change of the isobaric heat capacity of the solute moving from the mobile to the
stationary phase, in J mol⁻¹ K⁻¹.
* `φ₀`: Dimensionless film thickness (φ ≈ df/d) of the column for which the
thermodynamic parameters (Tchar, θchar, ΔCp) were estimated.
* `Cag`: diffusivity constant of solute `a` in gas `g`.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``H(x,t) = 2 \\frac{D_M}{u_M} + \\frac{d^2}{96}\\left(6 μ^2-16 μ +11
\\right) \\frac{u_M}{D_M} + \\frac{2}{3} d_f^2 μ(1-μ) \\frac{u_M}{D_S}``

with ``D_M`` the diffusion coefficient of the solute in the mobile phase,
``D_S`` the diffusion coefficient of the solute in the stationary phase,
``u_M`` the velocity of the mobile phase and μ the mobility of the solute.

``D_S`` is correlated to ``D_M`` by: 

``D_S = \\frac{D_M}{10000}``

**TODO**: alternative correlations?

``u_M`` is realated to the residency of the mobile phase ``r_M``:

``u_M = \\frac{1}{r_M}``

μ is correlated to the retention factor ``k``:

``μ = \\frac{1}{1 + k}``

See also: [`diffusion_mobile`](@ref), [`mobile_phase_residency`](@ref), [`retention_factor`](@ref)
"""
function plate_height(x, t, T_itp, Fpin_itp, pout_itp, L, d, df, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    id = d(x)# - 2.0*df(x)
    uM = 1/mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)
    μ = 1/(1 + retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀, k_th=k_th))
    DM = diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=ng, vis=vis, control=control)
    DS = DM/10000
    H1 = 2*DM/uM
    H2 = id^2/96*(6*μ^2-16*μ+11)*uM/DM
    H3 = 2/3*df(x)^2*μ*(1-μ)*uM/DS
    H = H1 + H2 + H3
    return H
end

function plate_height(x, t, T_itp, Fpin_itp, pout_itp, L, d::Number, df::Number, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    id = d# - 2.0*df
    uM = 1/mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)
    μ = 1/(1 + retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀, k_th=k_th))
    DM = diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=ng, vis=vis, control=control)
    DS = DM/10000
    H1 = 2*DM/uM
    H2 = id^2/96*(6*μ^2-16*μ+11)*uM/DM
    H3 = 2/3*df^2*μ*(1-μ)*uM/DS
    H = H1 + H2 + H3
    return H
end

function plate_height(x, t, T_itp, Fpin_itp, pout_itp, L, d, df::Number, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    id = d(x)# - 2.0*df
    uM = 1/mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)
    μ = 1/(1 + retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀, k_th=k_th))
    DM = diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=ng, vis=vis, control=control)
    DS = DM/10000
    H1 = 2*DM/uM
    H2 = id^2/96*(6*μ^2-16*μ+11)*uM/DM
    H3 = 2/3*df^2*μ*(1-μ)*uM/DS
    H = H1 + H2 + H3
    return H
end

function plate_height(x, t, T_itp, Fpin_itp, pout_itp, L, d::Number, df, gas, Tchar, θchar, ΔCp, φ₀, Cag; ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    id = d# - 2.0*df(x)
    uM = 1/mobile_phase_residency(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)
    μ = 1/(1 + retention_factor(x, t, T_itp, d, df, Tchar, θchar, ΔCp, φ₀, k_th=k_th))
    DM = diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=ng, vis=vis, control=control)
    DS = DM/10000
    H1 = 2*DM/uM
    H2 = id^2/96*(6*μ^2-16*μ+11)*uM/DM
    H3 = 2/3*df(x)^2*μ*(1-μ)*uM/DS
    H = H1 + H2 + H3
    return H
end

function plate_height(x, t, col::Column, prog::Program, sub::Substance, opt::Options)
    H = plate_height(x, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.df, col.gas, sub.Tchar, sub.θchar, sub.ΔCp, sub.φ₀, sub.Cag; ng=opt.ng, vis=opt.vis, control=opt.control, k_th=opt.k_th)
    return H
end

"""
    diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=false, vis="Blumberg", control="Pressure")

Calculate the diffusion coefficient of the solute in the mobile phase at
position `x` at time `t`.

# Arguments
* `x`: Position along the GC column, in m.
* `t`: Time in s.
* `T_itp`: Interpolated (linear) temperature `T(x,t)`.
* `Fpin_itp`: Interpolated (linear) Flow F(t) resp. inlet pressure `pin(t)`.
* `pout_itp`: Interpolated (linear) outlet pressure `pout(t)`.
* `L`: Length of the capillary measured in m (meter).
* `d`: Diameter of the GC column, in m.
* `gas`: Name of the mobile phase gas.
* `Cag`: diffusivity constant of solute `a` in gas `g`.
* `ng`: Option to calculate the simulation without a gradient (`ng = true`)
or with a gradient (`ng = false`).
* `vis`: used model for viscosity "Blumberg" or "HP".
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program

``D_M(x,t) = C_{ag} \\frac{T(x,t)^{1.75}}{p(x,t)}``
"""
function diffusion_mobile(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas, Cag; ng=false, vis="Blumberg", control="Pressure")
    DM = T_itp(x, t)^1.75/pressure(x, t, T_itp, Fpin_itp, pout_itp, L, d, gas; ng=ng, vis=vis, control=control)*Cag
    return DM
end

function diffusion_mobile(x, t, col::Column, prog::Program, sub::Substance, opt::Options)
    DM = diffusion_mobile(x, t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, col.L, col.d, col.gas, sub.Cag; ng=opt.ng, vis=opt.vis, control=opt.control)
    return DM
end

"""
    diffusivity(M, Cn, Hn, On, Nn, Rn, gas)

Calculate the diffusivity constant `Cag` of solute `a` in gas `g` using the
emperical Fuller-Schettler-Giddings model [1].

[1] Fuller, Edward N.; Ensley, Keith; Giddings, J. Calvin, Diffusion of
Halogenated Hydrocarbons in Helium. The Effect of Structure on Collision
Cross Sections, The Journal of Physical Chemistry, Volume 73, Issue 11,
1969, 3679–3685

# Arguments
* `M`: Molar mass of the solute.
* `Cn`: Number of carbon atoms of the solute.
* `Hn`: Number of hydrogen atoms of the solute.
* `On`: Number of oxygen atoms of the solute.
* `Nn`: Number of nitrogen atoms of the solute.
* `Rn`: Number of closed rings of the structure of the solute.
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.
"""
function diffusivity(M, Cn, Hn, On, Nn, Rn, gas)
    if gas=="H2"
        Vg = 6.12
        Mg = 2.02
    elseif gas=="He"
        Vg = 2.67
        Mg = 4
    elseif gas=="N2"
        Vg = 18.5
        Mg = 28.01
    elseif gas=="Ar"
        Vg = 16.2
        Mg = 39.95
    else
        error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
    end
    Va = 15.9*Cn + 2.31*Hn + 6.11*On + 4.54*Nn - 18.3*Rn
    Cag = pn*sqrt(1/M+1/Mg)/(Vg^(1/3)+Va^(1/3))^2*1e-7 # pn m²/s ('Cag at normal pressure')
    return Cag
end

"""
    diffusivity(id, gas)

Calculate the diffusivity constant `Cag` of solute `id` in gas `gas` using the
emperical Fuller-Schettler-Giddings model [1], using the identifier of the solute `id`.

[1] Fuller, Edward N.; Ensley, Keith; Giddings, J. Calvin, Diffusion of
Halogenated Hydrocarbons in Helium. The Effect of Structure on Collision
Cross Sections, The Journal of Physical Chemistry, Volume 73, Issue 11,
1969, 3679–3685

# Arguments
* `id`: Named tupel identifying the solute. (Name, CAS, formula, MW, smiles)
* `gas`: The name of the mobile phase. Allowed values: He, H2 or N2.
"""
function diffusivity(id, gas)
    if gas=="H2"
        Vg = 6.12
        Mg = 2.02
    elseif gas=="He"
        Vg = 2.67
        Mg = 4
    elseif gas=="N2"
        Vg = 18.5
        Mg = 28.01
    elseif gas=="Ar"
        Vg = 16.2
        Mg = 39.95
    else
        error("Unknown selection of gas. Choose one of these: He, H2, N2 or Ar.")
    end
    formula = formula_to_dict(id.formula)
    Rn = ring_number(id.smiles)
    Va = molecular_diffusion_volume(formula, Rn) # in cm³
    Cag = pn*sqrt(1/id.MW+1/Mg)/(Vg^(1/3)+Va^(1/3))^2*1e-7 # pn m²/s ('Cag at normal pressure')
    return Cag
end

"""
	formula_to_dict(formula)

Translate the formula string of a chemical substance into a dictionary, where the elements contained in the substance are the keys and the number of atoms are the values.

# Example 
```
julia> formula_to_dict("C14H20O")
Dict{String, Int64}("C" => 14, "H" => 20, "O" => 1)
```
"""
function formula_to_dict(formula)
	if ismissing(formula)
		formula_dict = missing
	else
		split_formula = eachmatch(r"[A-Z][a-z]*\d*", formula)
		elements = Array{Tuple{String, Int64}}(undef, length(collect(split_formula)))
		for i=1:length(collect(split_formula))
			formula_parts = collect(split_formula)[i].match
			element_string = match(r"[A-Za-z]+", formula_parts).match
			element_digit = match(r"[0-9]+", formula_parts)
			if isnothing(element_digit)
				elements[i] = (element_string, 1)
			else
				elements[i] = (element_string, parse(Int64, element_digit.match))
			end
		end
		formula_dict = Dict(elements)
	end
	return formula_dict
end

"""
	ring_number(smiles)

Extract the number of rings of a substance defined by its SMILES. The highest digit contained in the SMILES is returned as the number of rings. Only single digit ring numbers are recognized.
"""
function ring_number(smiles)
	if ismissing(smiles)
		rn = missing
	else
		allmatch = eachmatch(r"[0-9]", smiles)
		if isempty(allmatch)
			rn = 0
		else
			rn = maximum(parse.(Int64, [ match.match for match in allmatch]))
			# additional check, all integers from 1 to rn should be in allmatch, e.g. for rn=4 the integers 1, 2, 3, 4 should be in allmatch
		end
	end
	return rn
end

"""
    molecular_diffusion_volume(formula, Rn)

Calculate the molecular diffusion volume from the formula and number of rings of the solute as the sum of the atomic diffusion volumes according to [1].
    
[1] Fuller, Edward N.; Ensley, Keith; Giddings, J. Calvin, Diffusion of
Halogenated Hydrocarbons in Helium. The Effect of Structure on Collision
Cross Sections, The Journal of Physical Chemistry, Volume 73, Issue 11,
1969, 3679–3685

# Arguments
* `formula`: formula of the solute as a dictionary.
* `Rn`: number of rings of the solute.
"""
function molecular_diffusion_volume(formula, Rn) 
    if haskey(formula, "C")
        C = formula["C"]
    else
        C = 0
    end
    if haskey(formula, "H")
        H = formula["H"]
    else
        H = 0
    end
    if haskey(formula, "O")
        O = formula["O"]
    else
        O = 0
    end
    if haskey(formula, "N")
        N = formula["N"]
    else
        N = 0
    end
    if haskey(formula, "S")
        S = formula["S"]
    else
        S = 0
    end
    if haskey(formula, "F")
        F = formula["F"]
    else
        F = 0
    end
    if haskey(formula, "Cl")
        Cl = formula["Cl"]
    else
        Cl = 0
    end
    if haskey(formula, "Br")
        Br = formula["Br"]
    else
        Br = 0
    end
    if haskey(formula, "I")
        I = formula["I"]
    else
        I = 0
    end
    Va = C*15.9 + H*2.31 + O*6.11 + N*4.54 + S*22.9 + F*14.7 + Cl*21.0 + Br*21.9 + I*29.8 - Rn*18.3
    return Va
end

#---End-Functions-of-the-physical-model---