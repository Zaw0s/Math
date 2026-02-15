using Random, Distributions

function Categ(vec::Vector{<:Float64})
    k = length(vec)
    index = 1
    v = []
    for i in 1:k
        U = rand(Uniform(0,1))
        push!(v,-log(-log(U))+log(vec[i]))
    end
    for i in 2:k
        if v[index] < v[i]
            index = i
        end
    end
    return index
end

function evaluateNormal(x,mean,variance)
    return exp(-(x-mean)*(x-mean)/(2*variance))/sqrt(2*pi*variance)
end
function weightsForCateg(index,Y::Vector{<:Float64},θ::Vector{Vector{Float64}},w::Vector{<:Float64})
    m = length(w)
    weights = zeros(m)
    c = 0.0
    for i in 1: m
        c += evaluateNormal(Y[index],θ[i][1],1/θ[i][2])*w[i]
    end
    for i in 1:m
        weights[i] = evaluateNormal(Y[index],θ[i][1],1/θ[i][2])*w[i]
    end
    if c == 0.0
        return fill(1.0/m, m)
    end 
    
    return weights ./ c
end

function GibbsSamplerChain(Y::Vector{<:Float64},B,T,m,α::Vector{<:Float64},μ,λ,a,b)
    n = length(Y)
    d = zeros(Int, n)                       
    θ = [ [0.0, 1.0] for _ in 1:m ]         
    w = fill(1.0/m, m)                      

    for i in 1:n
        d[i] = rand(1:m)
    end

    Chain = []
    for i in 1:(B+T)
        for j in 1:n
            index = Categ(weightsForCateg(j,Y,θ,w))
            d[j] = index
        end
        for j in 1:m
            X = 0.0
            S = 0.0
            N = 0
            for k in 1:n
                if d[k] == j
                    X += Y[k]
                    N += 1
                end
            end
            if N > 0
                X/=N
                for k in 1:n
                    if d[k] == j
                        S += (Y[k]-X)^2
                    end
                end
                S/=N
            end

            μ1 = 0.0
            λ1 = 0.0
            α1 = 0.0
            β1 = 0.0

            if N > 0
                μ1 = (λ*μ+N*X)/(λ+N)
                λ1 = λ+N
                α1 = a+N/2
                β1 = b + (N*S+(λ*N*(X-μ)^2)/(λ+N))/2
            else
                μ1 = μ
                λ1 = λ
                α1 = a
                β1 = b 
            end
            #Creo que la definicion que nos dan es (shape,rate)
            τ = rand(Gamma(α1,1.0/β1))
            #julia usa desviacion estandar en lugar de varianza
            μ_τ = rand(Normal(μ1, sqrt(1/(λ1*τ))))

            θ[j][1] = μ_τ
            θ[j][2] = τ
        end

        simulationW = []
        for j in 1:m
            n_j = 0
            for k in 1:n
                if d[k] == j
                    n_j +=1
                end
            end
            #(shape,scale)
            push!(simulationW,rand(Gamma(α[j]+n_j,1)))
        end
        sumW = sum(simulationW)

        for j in 1:m
            w[j] = simulationW[j]/sumW
        end
        
        if i > B
            push!(Chain,[copy(d),deepcopy(θ),copy(w)])
        end
        
    end

    return Chain

end

data = parse.(Float64, readlines("galaxy.txt"))

MN = sum(data)/length(data)
C = GibbsSamplerChain(data,2000,10000,7,[1/7,1/7,1/7,1/7,1/7,1/7,1/7],MN,1/100,0.5,0.5)

function densityEstimator(x,MC)
    T = length(MC)
    m = length(MC[2])
    result = 0.0
    for i in 1:T
        for j in 1:m
            result += MC[i][3][j]*evaluateNormal(x,MC[i][2][1],1/MC[i][2][2])
        end
    end
    result /= T
    return result
end
