using LinearAlgebra

export fit_curve, normal, lognormal, exponential, powerlaw, linear, logistic, invlogistic

mutable struct normal
    x::Vector{<:Real}
    y::Vector{<:Real}
end
X(data::normal) =  hcat(ones(eltype(data.x), length(data.x)), data.x, (data.x).^2)
Y(data::normal) = log.(data.y)
invY(data::normal) = exp.(data.y)


mutable struct lognormal
    x::Vector{<:Real}
    y::Vector{<:Real}
end
X(data::lognormal) =  hcat(ones(eltype(data.x), length(data.x)), log.(data.x), (log.(data.x)).^2)
Y(data::lognormal) = log.(data.y)
invY(data::lognormal) = exp.(data.y)


mutable struct exponential
    x::Vector{<:Real}
    y::Vector{<:Real}
end
X(data::exponential) =  hcat(ones(eltype(data.x), length(data.x)), data.x)
Y(data::exponential) = log.(data.y)
invY(data::exponential) = exp.(data.y)


mutable struct powerlaw
    x::Vector{<:Real}
    y::Vector{<:Real}
end
X(data::powerlaw) =  hcat(ones(eltype(data.x), length(data.x)), log.(data.x))
Y(data::powerlaw) = log.(data.y)
invY(data::powerlaw) = exp.(data.y)


mutable struct linear
    x::Vector{<:Real}
    y::Vector{<:Real}
end
X(data::linear) =  hcat(ones(eltype(data.x), length(data.x)), data.x)
Y(data::linear) = data.y
invY(data::linear) = data.y


mutable struct logistic
    x::Vector{<:Real}
    y::Vector{<:Real}
end
X(data::logistic) =  hcat(ones(eltype(data.x), length(data.x)), data.x)
Y(data::logistic) = log.(inv.(data.y) .- 1.0)
invY(data::logistic) = inv.(exp.(data.y) .+ 1.0)


mutable struct invlogistic
    x::Vector{<:Real}
    y::Vector{<:Real}
end
X(data::invlogistic) =  hcat(ones(eltype(data.x), length(data.x)), data.x)
Y(data::invlogistic) = log.(1.0 .- data.y)
invY(data::invlogistic) = 1.0 .- exp.(data.y)


function smoothen(x::Vector{<:Real})
    x_min, x_max = extrema(x)
    return collect(range(x_min, stop = x_max, length = 5*length(x)))

end


function fit_curve(T::DataType, x::Vector{<:Real}, y::Vector{<:Real})
    data = T(x, y)
    
    fit_params = inv(adjoint(X(data))*X(data))*(adjoint(X(data))*Y(data))
    
    x_smooth = smoothen(x)
    data.x = x_smooth
    data.y = X(data)*fit_params
    y_smooth = invY(data)

    return x_smooth, y_smooth
end