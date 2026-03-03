# Fitting phi12 and phi13 together --> p1, p2, p3, K, C0, C1, C2, C3


using ADerrors, Optim, LaTeXStrings

include("../src/reader.jl")
include("../util/plot_utils.jl")
include("../util/err_utils.jl")

beta, phi2, phi12, phi13, phi4, fpik, a2t0, phik = Reader.retrieve_data()

########## Model ##########

function model12(phi2, a2t0, phi4, p)

    p1, p2, p3, K, C0, C1, C2, C3 = p

    phin = (4 .* phi4 .- 3 .* phi2) ./ 3

    L2 = phi2 .* log.(phi2)
    Ln = phin .* log.(phin)


    return phi2 .* (p1 .+ p2 .* phi2 .+ p3 .* K .* (L2 .- (1/3) .* Ln)) .+
           (a2t0 ./ 8) .* (C0 .+ C1 .* phi2)
end

function model13(phi2, a2t0, phi4, p)

    p1, p2, p3, K, C0, C1, C2, C3 = p

    phik = (2 .*phi4 .-phi2) ./2

    phin = (4 .* phi4 .- 3 .* phi2) ./ 3

    Ln = phin .* log.(phin)


    return (phik.*(p1 .+p2 .*phik .+ 2/3 .* p3.*K.*Ln) + a2t0 .*(C2 .+ C3 .*phi2) ./8)
end

function model(phi2, a2t0, phi4, p)
    return vcat(model12(phi2, a2t0, phi4, p), model13(phi2, a2t0, phi4, p))
end

function chi2_fit12(params)
    theory = model12(value.(phi2), value.(a2t0), value.(phi4), params)
    return sum( (value.(phi12)-theory).^2 ./ (ADerrors.err.(phi12)).^2 )
end

function chi2_fit13(params)
    theory = model13(value.(phi2), value.(a2t0), value.(phi4), params)
    return sum( (value.(phi13)-theory).^2 ./ (ADerrors.err.(phi13)).^2 )
end

function chi2_fit(params)
    return (chi2_fit12(params) + chi2_fit13(params))
end

p0 = [0.13, 0.0, 0.01, 1.0, 0.0, 0.0, 0.0, 0.0]
result = optimize(chi2_fit, p0, BFGS())
params = Optim.minimizer(result)
# [0.1424381908921298, -0.002019464229222623, 0.009840659846048896, 1.000098078763568, -0.0362720876284107, -0.3124277427543239, -0.4933508451266846, 0.3241243146078352]

### Propagating error to params ###

function chisq(p, d)
    theory = model(value.(phi2), value.(a2t0), value.(phi4), p)
    err = vcat(ADerrors.err.(phi12), ADerrors.err.(phi13))
    return sum((d .- theory).^2 ./ err.^2)
end

fit_data = vcat(phi12, phi13)
(fitp, cexp) = fit_error(chisq, params, fit_data)


ErrUtils.LoadErr(fitp)

### Plotting bands ###

#= function model_plot12(phi2, a2t0, phi4, p)

    p1, p2, p3, K, C0, C1, C2, C3 = p

    phin = (4 * phi4 - 3 * phi2) / 3

    L2 = phi2 * log(phi2)
    Ln = phin * log(phin)

    return phi2 * (p1 + p2 * phi2 + p3 * K * (L2 - (1/3) * Ln)) +
           (a2t0 / 8) * (C0 + C1 * phi2)
end

xaxis = range(minimum(value.(phi2)), maximum(value.(phi2)), length=300)
bs = [3.4, 3.46, 3.55, 3.7, 3.85]
as = [a2t0[1], a2t0[4], a2t0[6], a2t0[11], a2t0[15]] # Lattice spacings associated to each beta

for i in eachindex(bs)
    band = uwreal[]

    for x in xaxis
        val = model_plot12(x, as[i], phi4[1], fitp)
        push!(band, val)
    end

    y12 = band./xaxis
    ErrUtils.LoadErr(y12)

    PlotUtils.CreateBand(xaxis, y12, PlotUtils.colors[i])
end

# Continuum band

band = uwreal[]

for x in xaxis
        val = model_plot12(x, 0, phi4[1], fitp)
        push!(band, val)
    end

    y12 = band./xaxis
    ErrUtils.LoadErr(y12)

    PlotUtils.CreateBand(xaxis, y12, :gray, "Continuum")

# Scatter data

y = phi12./phi2
ErrUtils.LoadErr(y)
PlotUtils.ScatterByEnsemble(beta, phi2, y, L"\phi_2", L"\phi_{12}/\phi_2", latexstring("\\mathrm{Joint\\ fit\\ of}\\ \\phi_{12}\\ \\mathrm{and}\\ \\phi_{13}"), "phi12_phi13/phi12")
 =#


function model_plot13(phi2, a2t0, phi4, p)

    p1, p2, p3, K, C0, C1, C2, C3 = p

    phik = (2 *phi4 -phi2) /2

    phin = (4 * phi4 - 3 * phi2) / 3

    Ln = phin * log(phin)


    return (phik*(p1 +p2 *phik + 2/3 * p3*K*Ln) + a2t0 *(C2 + C3 *phi2) /8)
end

xaxis = range(minimum(value.(phi2)), maximum(value.(phi2)), length=300)
phik_plot = [(2*phi4[1]-x)/2 for x in xaxis]

bs = [3.4, 3.46, 3.55, 3.7, 3.85]
as = [a2t0[1], a2t0[4], a2t0[6], a2t0[11], a2t0[15]] # Lattice spacings associated to each beta

for i in eachindex(bs)
    band = uwreal[]

    for x in xaxis
        val = model_plot13(x, as[i], phi4[1], fitp)
        push!(band, val)
    end

    y13 = band./phik_plot
    ErrUtils.LoadErr(y13)

    PlotUtils.CreateBand(xaxis, y13, PlotUtils.colors[i])
end

# Continuum band

band = uwreal[]

for x in xaxis
        val = model_plot13(x, 0, phi4[1], fitp)
        push!(band, val)
    end

    y13 = band./phik_plot
    ErrUtils.LoadErr(y13)

    PlotUtils.CreateBand(xaxis, y13, :gray, "Continuum")

# Scatter data

y = phi13./phik
ErrUtils.LoadErr(y)
PlotUtils.ScatterByEnsemble(beta, phi2, y, L"\phi_2", L"\phi_{13}/\phi_K", latexstring("\\mathrm{Joint\\ fit\\ of}\\ \\phi_{12}\\ \\mathrm{and}\\ \\phi_{13}"), "phi12_phi13/phi13")

### Dumping details of params ###

ErrUtils.DumpErr(fitp, ["p1", "p2", "p3", "K", "C0", "C1", "C'0 (C2)", "C'1 (C3)"], "phi12_phi13")