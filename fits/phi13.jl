# Fitting phi13(phi2) --> p1, p2, p3, K, C0, C1

using ADerrors, Optim, LaTeXStrings

include("../src/reader.jl")
include("../util/plot_utils.jl")
include("../util/err_utils.jl")

beta, phi2, phi12, phi13, phi4, fpik, a2t0, phik = Reader.retrieve_data()

########## Model ##########

function model(phi2, a2t0, phi4, p)

    p1, p2, p3, K, C0, C1 = p

    phik = (2 .*phi4 .-phi2) ./2

    phin = (4 .* phi4 .- 3 .* phi2) ./ 3

    Ln = phin .* log.(phin)


    return (phik.*(p1 .+p2 .*phik .+ 2/3 .* p3.*K.*Ln) + a2t0 .*(C0 .+ C1 .*phi2) ./8)
end

function chi2_fit(params)
    theory = model(value.(phi2), value.(a2t0), value.(phi4), params)
    return sum( (value.(phi13)-theory).^2 ./ (ADerrors.err.(phi13)).^2 )
end

p0 = [0.13, 0.0, 0.01, 1.0, 0.0, 0.0]
result = optimize(chi2_fit, p0, BFGS())
params = Optim.minimizer(result)
# [0.0929011444908818, 0.05686999282293377, -0.028377699380523046, 0.9706256552110858, -0.579990774676675, 0.4755241779774559]

### Propagating error to params ###

chisq(p, d) = sum( (d .- model(value.(phi2), value.(a2t0), value.(phi4), p)) .^2 ./ (ADerrors.err.(phi13)).^2)

(fitp, cexp) = fit_error(chisq, params, phi13)

ErrUtils.LoadErr(fitp)

### Plotting bands ###

function model_plot(phi2, a2t0, phi4, p)

    p1, p2, p3, K, C0, C1 = p

    phik = (2 *phi4 -phi2) /2

    phin = (4 * phi4 - 3 * phi2) / 3

    Ln = phin * log(phin)


    return (phik*(p1 +p2 *phik + 2/3 * p3*K*Ln) + a2t0 *(C0 + C1 *phi2) /8)
end

xaxis = range(minimum(value.(phi2)), maximum(value.(phi2)), length=300)
phik_plot = [(2*phi4[1]-x)/2 for x in xaxis]

bs = [3.4, 3.46, 3.55, 3.7, 3.85]
as = [a2t0[1], a2t0[4], a2t0[6], a2t0[11], a2t0[15]] # Lattice spacings associated to each beta

for i in eachindex(bs)
    band = uwreal[]

    for x in xaxis
        val = model_plot(x, as[i], phi4[1], fitp)
        push!(band, val)
    end

    y13 = band./phik_plot
    ErrUtils.LoadErr(y13)

    PlotUtils.CreateBand(xaxis, y13, PlotUtils.colors[i])
end

# Continuum band

band = uwreal[]

for x in xaxis
        val = model_plot(x, 0, phi4[1], fitp)
        push!(band, val)
    end

    y13 = band./phik_plot
    ErrUtils.LoadErr(y13)

    PlotUtils.CreateBand(xaxis, y13, :gray, "Continuum")

# Scatter data

y = phi13./phik
ErrUtils.LoadErr(y)
PlotUtils.ScatterByEnsemble(beta, phi2, y, L"\phi_2", L"\phi_{13}/\phi_K", latexstring("\\mathrm{Single\\ fit\\ of}\\ \\phi_{13}"), "phi13/phi13")

### Dumping details of params ###

ErrUtils.DumpErr(fitp, ["p1", "p2", "p3", "K", "C'0", "C'1"], "phi13")