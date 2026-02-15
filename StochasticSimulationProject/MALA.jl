using Random, Distributions, Plots, DataStructures, SpecialFunctions,QuadGK
function S()
    return rand(Normal(0,1))
end
function U()
    return rand(Uniform(0,1))
end

function V(x)
    return 0.25*x^4-4.5*x^2
end

function ∇V(x)
    return x^3-9.0*x
end

function a(x,y,h)
    val = -V(y) + V(x) + ((y-x+h*∇V(x))^2-(x-y+h*∇V(y))^2)/(4.0*h)
    val = min(1.0,exp(val))
    return val
end

function ULA(n,h)
    X =  0.0
    first = Float64[0.0]
    last = Float64[]
    p = Float64[0.0] 
    for i in 1:n
        X = X-h*∇V(X)+sqrt(2.0*h)*S()
        push!(p,X)
        if i <= 1999
            push!(first,X)
        end
        if i >= n-1999
            push!(last,X)
        end
    end
    #x = 1:length(first)
    #plot(x,first, label = "primeras 2000")
    #plot!(x,last, label = "últimas 2000")
    plot(p)
    #histogram(p)
end

ULA(1000000,0.1)

function MALA(n,h)
    X =  0.0
    first = [0.0]
    last = []
    p = [0.0]
    for i in 1:n
        Y = X-h*∇V(X)+sqrt(2*h)*S()
        if U() < a(X,Y,h)
            X=Y
        end
        push!(p,X)
        if i <= 1999
            push!(first,X)
        end
        if i >= n-1999
            push!(last,X)
        end
    end
    #x = 1:length(first)
    #plot(x,first, label = "primeras 2000")
    #plot!(x,last, label = "últimas 2000")
    plot(p)
    #histogram(p)
end

MALA(1000900,0.1)

