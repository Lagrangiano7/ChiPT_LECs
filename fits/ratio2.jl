# Fitting phi12(phi2) --> p1, p2, p3, K, G0, G1

using ADerrors, Optim, LaTeXStrings

include("../src/reader.jl")
include("../util/plot_utils.jl")
include("../util/err_utils.jl")

beta, phi2, phi12, phi13, phi4, fpik, a2t0, phik = Reader.retrieve_data()

ratio2 = 4 .*phi13 ./(2 .*phi4 .- phi2) + phi12./phi2
ErrUtils.LoadErr(ratio2)

########## Model ##########

function model(phi2, a2t0, phi4, p)

    p1, p2, p3, K, G0, G1 = p

    phin = (4 .* phi4 .- 3 .* phi2) ./ 3

    L2 = phi2 .* log.(phi2)
    Ln = phin .* log.(phin)


    return 3 .* p1 .+ 2 .* p2 .* phi4 .+ p3 .* K.* (L2 .- Ln) .+ a2t0 .* (G0 .+ G1 .* phi2) ./8
end

function chi2_fit(params)
    theory = model(value.(phi2), value.(a2t0), value.(phi4), params)
    return sum( (value.(ratio2)-theory).^2 ./ (ADerrors.err.(ratio2)).^2 )
end

p0 = [0.13, 0.0, 0.01, 1.0, 0.0, 0.0]
result = optimize(chi2_fit, p0, BFGS())
params = Optim.minimizer(result)

### Propagating error to params ###

chisq(p, d) = sum( (d .- model(value.(phi2), value.(a2t0), value.(phi4), p)) .^2 ./ (ADerrors.err.(ratio2)).^2)

(fitp, cexp) = fit_error(chisq, params, ratio2)

ErrUtils.LoadErr(fitp)

### Plotting bands ###

function model_plot(phi2, a2t0, phi4, p)

    p1, p2, p3, K, G0, G1 = p

    phin = (4 * phi4 - 3 * phi2) / 3

    L2 = phi2 * log(phi2)
    Ln = phin * log(phin)


    return 3 *p1 + 2 *p2 *phi4 + p3*K*(L2 - Ln) + a2t0 *(G0 +G1*phi2) /8
end

xaxis = range(minimum(value.(phi2)), maximum(value.(phi2)), length=300)
bs = [3.4, 3.46, 3.55, 3.7, 3.85]
as = [a2t0[1], a2t0[4], a2t0[6], a2t0[11], a2t0[15]] # Lattice spacings associated to each beta

for i in eachindex(bs)
    band = uwreal[]

    for x in xaxis
        val = model_plot(x, as[i], phi4[1], fitp)
        push!(band, val)
    end

    y12 = band
    ErrUtils.LoadErr(y12)

    PlotUtils.CreateBand(xaxis, y12, PlotUtils.colors[i])
end

# Continuum band

band = uwreal[]

for x in xaxis
        val = model_plot(x, 0, phi4[1], fitp)
        push!(band, val)
    end

    y12 = band
    ErrUtils.LoadErr(y12)

    PlotUtils.CreateBand(xaxis, y12, :gray, "Continuum")

# Scatter data


PlotUtils.ScatterByEnsemble(beta, phi2, ratio2, L"\phi_2", L"\frac{4\phi_{13}}{2\phi_4-\phi_2}+\frac{\phi_{12}}{\phi_2}", "", "ratio2/ratio2", true)

### Dumping details of params ###

ErrUtils.DumpErr(fitp, ["p1", "p2", "p3", "K", "G0", "G1"], "ratio2")