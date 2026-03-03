module Reader

using ADerrors, BDIO, ALPHAio, LsqFit, Statistics, Optim, Plots
include("./const.jl")

######## READING FILE ########
fn = "C:\\Users\\germa\\Desktop\\Practicas\\src\\mq_alpha.bdio"
fb = BDIO_open(fn, "r")
obs = []
dobs = []
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    push!(dobs, d)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)

    dw = ALPHAdobs_read_next(fb, size=sz)
    push!(obs, dw)
    #dobs -> dict and info, extra
    #e.g. dobs[1]["extra"]["obs"] -> "am12"
    #e.g. dobs[1]["extra"]["ens_id"][2] -> "H400"
    #obs  -> data
    #e.g. obs[1][2] -> am12[H400]
end
BDIO_close!(fb)
m12 = [ntuple(i -> obs[1][i],16)...]
m13 = [ntuple(i -> obs[2][i],16)...]
beta = dobs[1]["extra"]["beta"]
phi13 = [ntuple(i -> m13[i]*sqrt(8*obs[5][i]),16)...]
phi12 = [ntuple(i -> m12[i]*sqrt(8*obs[5][i]),16)...]
phi2 = [ntuple(i -> obs[3][i],16)...]
phi4 = [ntuple(i -> obs[4][i],16)...]
a2t0 = [ntuple(i -> 1/obs[5][i],16)...]

fpik = [ntuple(i -> obs[6][i]/sqrt(obs[5][i]),16)...]
ens_list = [ntuple(i -> dobs[1]["extra"]["ens_id"][i],16)...]
phi4ref = phi4[1].mean
phi2symval = phi2[1].mean
phi2sym = phi2[1]

######## End of reading file ########

phik = (2 .*phi4 .-phi2) ./2

for i in eachindex(phi2)
    uwerr(phi12[i])
    uwerr(phi13[i])
    uwerr(phi2[i])
    uwerr(phi4[i])
    uwerr(a2t0[i])
    uwerr(fpik[i])

    uwerr(phik[i])
end

function retrieve_data()
    return (beta, phi2, phi12, phi13, phi4, fpik, a2t0, phik)
end



end