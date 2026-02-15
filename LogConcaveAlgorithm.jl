using Random, Distributions, Plots, DataStructures, SpecialFunctions,QuadGK
function evaluateNormal(x)
    mean = 0
    variance = 1
    return exp(-(x-mean)*(x-mean)/(2*variance))/sqrt(2*pi*variance)
end
function derMinusLogNormal(x)
    return x
end

function evaluateGamma(x,a)
    if x <= -(a-1)
        return 0
    end
    return  ((x+(a-1))^(a-1))*exp(-(x+(a-1)))/gamma(a)
end

function derMinusLogGamma(x,a)
    return -(a-1)/(x+(a-1))+1
end


function evaluateBeta(x,a,b)
    if x <= -(a-1)/(a+b-2) || x >= 1-(a-1)/(a+b-2)
        return 0 
    end
    return (((x+(a-1)/(a+b-2))^(a-1))*(1-(x+(a-1)/(a+b-2)))^(b-1))/beta(a,b)
end

function derMinusLogBeta(x,a,b)
    return -(a-1)/(x+(a-1)/(a+b-2))+(b-1)/(1-(x+(a-1)/(a+b-2)))
end

function E(x,type,α,Β)
    if type == 1
        return evaluateNormal(x)
    elseif type == 2
        return evaluateGamma(x,α)
    else
        return evaluateBeta(x,α,Β)
    end
end
function D(x,type,α,β)
    if type == 1
        return derMinusLogNormal(x)
    elseif type == 2
        return derMinusLogGamma(x,α)
    else
        return derMinusLogBeta(x,α,β)
    end
end
function firstAlgorithm(E,q,n,type,a,b)
    c = 2
    simulation = []
    ξ = E(0,type,a,b)#f(0)
    while length(simulation) < n
        U = rand()
        X = 0
        val = 0
        
        if U < q/c
            X = (log(U*c/q)-1)*q/ξ
            val = exp(1+X*ξ/q)*ξ
        elseif U < (1+q)/c
            X = (c*U-2*q)/ξ
            val = ξ
        else
            X = (1 - log((c-U*c)/(1-q)))*(1-q)/ξ
            val = exp(1-X*ξ/(1-q))*ξ
        end

        U1 = rand()
        if U1 <= E(X,type,a,b)/(val)
            push!(simulation,X)
        end
    end
    return simulation
end

histogram(firstAlgorithm(E,1/2,1000000,1,0,0))

histogram(firstAlgorithm(E,1-gamma(3,2)/gamma(3),1000000,2,3,0))
histogram(firstAlgorithm(E,1-gamma(4,3)/gamma(4),1000000,2,4,0))
histogram(firstAlgorithm(E,1-gamma(5,4)/gamma(5),1000000,2,5,0))
histogram(firstAlgorithm(E,1-gamma(6,5)/gamma(6),1000000,2,6,0))
histogram(firstAlgorithm(E,1-gamma(7,6)/gamma(7),1000000,2,7,0))



histogram(firstAlgorithm(E,beta_inc(2, 2,1/2)[1],1000000,3,2,2))
histogram(firstAlgorithm(E,beta_inc(1.5, 1.5,1/2)[1],1000000,3,1.5,1.5))
histogram(firstAlgorithm(E,beta_inc(3, 3,1/2)[1],1000000,3,3,3))
histogram(firstAlgorithm(E,beta_inc(2.2, 2,1.2/(2.2))[1],1000000,3,2.2,2))
histogram(firstAlgorithm(E,beta_inc(1.5, 2.2,0.5/(1.7))[1],1000000,3,1.5,2.2))

#derivate of -log(f)
function secondAlgorithm(E,D,n,typeOfFunction,α,β)
    posDerivate = SortedDict{Float64,Float64}()

    insert!(posDerivate,-0.2,D(-0.2,typeOfFunction,α,β))
    insert!(posDerivate,0,D(0,typeOfFunction,α,β))
    insert!(posDerivate,0.2,D(0.2,typeOfFunction,α,β))
        
    simulation = []
    
    while length(simulation) < n

        intersection = []
        
        flag = 0
        x = 0
        m = 0

        for (pos,der) in posDerivate
            if flag == 0
                flag = 1
            else
                b1 = -log(E(x,typeOfFunction,α,β))-m*x
                b2 = -log(E(pos,typeOfFunction,α,β))-der*pos
                push!(intersection,(b2-b1)/(m-der))
            end
            x = pos
            m = der
        end
        index = 1
        integral = 0.0
        for (pos,der) in posDerivate
            b = 0
            if der == 0.0                    
                b = -log(E(pos,typeOfFunction,α,β))-der*pos
                integral += (intersection[index]-intersection[index-1])*exp(-b)
            else
                if index > 1
                    b = -log(E(pos,typeOfFunction,α,β))-der*pos
                    integral += exp(-der*intersection[index-1]-b)/der 
                end
                if index <= length(intersection)
                    b = -log(E(pos,typeOfFunction,α,β))-der*pos
                    integral -= exp(-der*intersection[index]-b)/der        
                end
            end
            index += 1
        end
        index = 1
        acc = 0.0
        U = rand(Uniform(0,integral))
        X = 0
        h = 0

        for (pos,der) in posDerivate
            intervalIntegral = 0
            if der == 0.0                    
                b = -log(E(pos,typeOfFunction,α,β))-der*pos
                intervalIntegral += (intersection[index]-intersection[index-1])*exp(-b)
                if U < acc +intervalIntegral
                    X = ((U-acc)/exp(-b))+intersection[index-1]
                    h = exp(-der*X-b)
                    break
                end
            else
                if index == 1

                    b = -log(E(pos,typeOfFunction,α,β))-der*pos
                    intervalIntegral -= exp(-der*intersection[index]-b)/der

                    if U < acc +intervalIntegral
                        X = (b +log(-der*U))/(-der)
                        h = exp(-der*X-b)
                        break
                    end
                elseif index <= length(intersection)

                    b = -log(E(pos,typeOfFunction,α,β))-der*pos
                    intervalIntegral -= exp(-der*intersection[index]-b)/der
                    intervalIntegral += exp(-der*intersection[index-1]-b)/der

                    if U < acc +intervalIntegral
                        X = (b +log((-der*U+der*acc+exp(-der*intersection[index-1]-b))))/(-der)
                        h = exp(-der*X-b)                   
                        break
                    end
                else

                    b = -log(E(pos,typeOfFunction,α,β))-der*pos
                    intervalIntegral += exp(-der*intersection[index-1]-b)/der

                    X = (b +log((-der*U+der*acc+exp(-der*intersection[index-1]-b))))/(-der)
                    h = exp(-der*X-b)                   
                end
            end
            acc+=intervalIntegral
            index += 1
        end

        U1 = rand(Uniform(0,1))
        
        if U1 < E(X,typeOfFunction,α,β)/h
            push!(simulation,X)
        else
            if abs(X) < 10 && E(X,typeOfFunction,α,β) > 0.05
                insert!(posDerivate,X,D(X,typeOfFunction,α,β))
            end
        end

    end
    histogram(simulation)
    
end

secondAlgorithm(E,D,100000,1,0,0)

secondAlgorithm(E,D,100000,2,3,0)
secondAlgorithm(E,D,100000,2,4,0)
secondAlgorithm(E,D,100000,2,5,0)
secondAlgorithm(E,D,100000,2,6,0)
secondAlgorithm(E,D,100000,2,7,0)

secondAlgorithm(E,D,100000,3,2,2)
secondAlgorithm(E,D,100000,3,1.5,1.5)
secondAlgorithm(E,D,100000,3,3,3)
secondAlgorithm(E,D,100000,3,2.2,2)
secondAlgorithm(E,D,100000,3,1.5,2.2)

