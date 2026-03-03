module PlotUtils

using Plots, ADerrors

colors = (:red, :blue, :green, :black, :orange)

function ScatterByEnsemble(beta, xdata, ydata, xlabel="", ylabel="", title="", filename="default", enhanceLegend=false)
    xdata_val = value.(xdata)
    ydata_val = value.(ydata)
    ydata_err = ADerrors.err.(ydata)

    color_ind = 1
    cur_beta = beta[1]

    # Physical value of phi2
    vline!([minimum(xdata_val)], linestyle=:dash, color=:black, xlims=[0.0, maximum(xdata_val)*1.1], label="", grid=false)

    if enhanceLegend
        ylims = extrema(ydata_val)
        display(scatter!((xdata_val[1], ydata_val[1]), xlabel=xlabel, ylabel=ylabel, title=title, color=colors[color_ind], yerr=ydata_err[1], label="β=$cur_beta", ylim=(ylims[1]*0.9, :auto), legend=:bottomright))

    else
        display(scatter!((xdata_val[1], ydata_val[1]), xlabel=xlabel, ylabel=ylabel, title=title, color=colors[color_ind], yerr=ydata_err[1], label="β=$cur_beta"))
    end

    for i in eachindex(beta)
        # Update beta and its corresponding colour code if proceeding
        if i==1
            continue
        end
        if beta[i] != cur_beta
            color_ind += 1
            cur_beta = beta[i]
            display(scatter!((xdata_val[i], ydata_val[i]), color=colors[color_ind], yerr=ydata_err[i], label="β=$cur_beta"))
        else
            display(scatter!((xdata_val[i], ydata_val[i]), color=colors[color_ind], yerr=ydata_err[i], label=""))
        end
        
    end

    if filename!="default"
        savefig("../img/$filename.png")
    else
        println("Press any key to end")
        readline()
    end
end

function CreateBand(xaxis, yaxis, c, label="")
    plot!(xaxis, value.(yaxis),
     ribbon=ADerrors.err.(yaxis),
     fillalpha=0.25,
     color=c,
     lw=2,
     label=label)
end

end