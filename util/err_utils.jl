module ErrUtils

using ADerrors

function LoadErr(x)
    for i in eachindex(x)
        uwerr(x[i])
    end
end

function DumpErr(params, names, filename)
    open("../fits/details/$filename.txt", "w") do file
        redirect_stdout(file) do
            for i in eachindex(params)
                println(names[i])
                details(params[i])
                println()
            end
        end
    end
end


### End of Module ###
end