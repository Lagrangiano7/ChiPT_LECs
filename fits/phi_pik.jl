# Fitting phi_pik(phi2) --> A, B, C0, C1

using ADerrors, Optim, LaTeXStrings

include("../src/reader.jl")
include("../util/plot_utils.jl")
include("../util/err_utils.jl")

beta, phi2, phi12, phi13, phi4, fpik, a2t0, phik = Reader.retrieve_data()

phi_pik = sqrt(8)*fpik./sqrt.(a2t0)
ErrUtils.LoadErr(phi_pik)

########## Model ##########

L(x) = x*log(x)

function model(phi2, a2t0, phi4, p)
    A, B, C0, C1 = p

    return A ./ (4 .*pi) .* (1 .- 7/6 .*L.(phi2./(A.^2)) - 4/3 .* L.((phi4 .- 1/2 .*phi2)./(A.^2)) - 1/2 .* L.((4/3 .*phi4 .- phi2)./(A.^2)) + B./(A.^2) .*phi4) + a2t0.*(C0 .+ C1.*phi2) ./ 8
end

function chi2_fit(params)
    theory = model(value.(phi2), value.(a2t0), value.(phi4), params)
    return sum( (value.(phi_pik)-theory).^2 ./ (ADerrors.err.(phi_pik)).^2 )
end

p0 = [4.8, 1.0, 0.0, 0.0]
result = optimize(chi2_fit, p0, BFGS())
params = Optim.minimizer(result)

### Propagating error to params ###

chisq(p, d) = sum( (d .- model(value.(phi2), value.(a2t0), value.(phi4), p)) .^2 ./ (ADerrors.err.(phi_pik)).^2)

(fitp, cexp) = fit_error(chisq, params, phi_pik)

ErrUtils.LoadErr(fitp)

### Plotting bands ###

function model_plot(phi2, a2t0, phi4, p)
    A, B, C0, C1 = p

    return A / (4*pi) * (1 - 7/6 *L(phi2/(A^2)) - 4/3* L((phi4 - 1/2 *phi2)/(A^2)) - 1/2* L((4/3 *phi4 - phi2)/(A^2)) + B/(A^2) *phi4) + a2t0*(C0 + C1*phi2) / 8
end

xaxis = range(minimum(value.(phi2)*1e-3), maximum(value.(phi2)), length=300)

bs = [3.4, 3.46, 3.55, 3.7, 3.85]
as = [a2t0[1], a2t0[4], a2t0[6], a2t0[11], a2t0[15]] # Lattice spacings associated to each beta

for i in eachindex(bs)
    band = uwreal[]

    for x in xaxis
        val = model_plot(x, as[i], phi4[1], fitp)
        push!(band, val)
    end

    ypik = band
    ErrUtils.LoadErr(ypik)

    PlotUtils.CreateBand(xaxis, ypik, PlotUtils.colors[i])
end

# Continuum band

band = uwreal[]

for x in xaxis
    val = model_plot(x, 0, phi4[1], fitp)
    push!(band, val)
end

ypik = band
ErrUtils.LoadErr(ypik)

PlotUtils.CreateBand(xaxis, ypik, :gray, "Continuum")

# Scatter data


PlotUtils.ScatterByEnsemble(beta, phi2, phi_pik, L"\phi_2", L"\phi_{\pi K}", L"Single fit of $\phi_{\pi K}$", "phi_pik/phi_pik")

### Dumping details of params ###

ErrUtils.DumpErr(fitp, ["A", "B", "C0", "C1"], "phi_pik")