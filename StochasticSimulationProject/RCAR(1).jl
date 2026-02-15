using Random, Distributions, Plots, DataStructures, SpecialFunctions,QuadGK

function S()
    return rand(Normal(0,1))
end
function U()
    return rand(Uniform(0,1))
end
function MCMC(S,U,x0::Float64,steps)
    μ = x0
    X = x0
    simulation = [x0]
    estim = [x0]
    for i in 1:steps
        X = X*U()+S()
        μ = μ*(i/(i+1))+X/(i+1)
        push!(estim,last(estim)+X)
        push!(simulation,X)
    end
    xtime = []
    for i in 0:steps
        estim[i+1]/=(i+1)
        push!(xtime,i)
    end
    p = plot(xtime,estim)
    plot!(p,xtime,simulation)
    display(p)
end

MCMC(S,U,100.0, 100)

function MLMC(S,U,x0::Float64,e::Float64)
    C = abs(x0)
    if C == 0
        println("ECM es 0")
    else
        print(ceil((log(C)-log(e/sqrt(2)))/(-log(0.5))))
        steps =  Int64(ceil((log(C)-log(e/sqrt(2)))/(-log(0.5))))
    end
    V = [0.0]
    N = [1]
    f = 1.0
    for i in 1:steps
        calc = 3*(abs(x0)*f/2+abs(x0))
        push!(V,calc)
        f*=0.5
    end
    k = 0.0
    for i in 1:steps+1
        k += sqrt(V[i]/i)
    end
    k *= 2.0/(e^2)

    for i in 1:steps
        push!(N,ceil(k*sqrt(V[i+1]/(i+1))))
    end

    simulation =[]
    for i in 0:steps
        sum = 0.0
        for j in 0: i
            estim = 0.0
            for samp in 1: N[j+1]
                u = 1.0
                for r1 in 1:(j-1)
                    u *= U()
                end
                if j == 0
                    estim += x0
                else
                    estim += u*(x0*U()-x0+S())
                end
            end
            if N[j+1] > 0
                estim/=N[j+1]
            end
            sum += estim
        end
        push!(simulation,sum)
    end

    xtime = []
    for i in 0:steps
        push!(xtime,i)
    end
    plot(xtime,simulation)
end

MLMC(S,U,-50.0,0.1)
