using StatsBase
using Statistics
using LinearAlgebra
using MultivariateStats
using FileIO, JLD

export load_jld_data, col_normalize, reverse_matrix, smooth_line, auto_reshape, get_cdist, model_sigmoid, model_poly, model_rc, get_dist, print_progress

function load_jld_data(filename::String, varname::String)
    isdir("data-files") ? fname = string(pwd(), "/data-files/", filename) : error("data-files not in working directory")
    
    if isfile(fname)
        return load(File{format"JLD"}(fname), varname)
    else
        error("file not in data-files")
    end
end

function col_normalize(x::AbstractArray{<:Real}; showNorms::Bool = false)
    if isa(x, Union{Matrix, LinearAlgebra.Adjoint, SubArray{eltype(x), 2}})
        m, n = size(x)
    elseif isa(x, Union{Vector, SubArray{eltype(x), 1}})
        m = length(x)
        n = 1
    end

    y = Array{Float64, 2}(undef, m, n)
    col_norms = Array{Float64, 1}(undef, n)

    for i in 1:n
        col = @view x[:, i]
        @inbounds col_norms[i] = norm(col)
        @inbounds col_norms[i] != 0.0 ? y[:, i] .= col./col_norms[i] : y[:, i] .= 0.0
    end

    if showNorms 
        return y, col_norms
    else
        return y
    end
    
end;

function col_normalize!(dest::AbstractArray{Float64}, x::AbstractArray{<:Real}; showNorms::Bool = false)
    if isa(x, Union{Matrix, LinearAlgebra.Adjoint, SubArray{eltype(x), 2}})
        n = size(x, 2)
    elseif isa(x, Union{Vector, SubArray{eltype(x), 1}})
        n = 1
    end

    col_norms = Array{Float64, 1}(undef, n)

    for i in 1:n
        col = @view x[:, i]
        @inbounds col_norms[i] = norm(col)
        @inbounds col_norms[i] != 0.0 ? dest[:, i] .= col./col_norms[i] : dest[:, i] .= 0.0
    end

    if showNorms
        return col_norms
    end
    
end;

function reverse_matrix(matrix::AbstractArray{<:Real, 2})
    (nrows, ncols) = size(matrix)
    x = zeros(size(matrix))
    for i = 1:nrows
        x[i,:] = @view matrix[nrows+1-i,:]
    end
    
    return x
end;


function smooth_line(data::Vector{<:Real}, window_size::Int64)
    y_smooth = zeros(length(data))
    
    for index = 1:length(data)
        if(index < window_size)
            y_smooth[index] = mean(data[1:index])
        end
        if(index >= window_size)
            y_smooth[index] = mean(data[index - window_size + 1:index])
        end
    end
    
    return y_smooth
end;

function auto_reshape(x::AbstractArray{<:Real, 1})
    array_length = length(x)
    nrows = 0
    ncols = 0
    
    for delta = 0:array_length-1
        nrows = 0.5*(-delta + sqrt(delta^2 + 4*array_length))
        
        if rem(nrows, 1) == 0
            ncols = array_length/nrows
            break
        end
    end
    
    return reshape(x, convert(Int64, nrows), convert(Int64, ncols))
end;

function get_cdist(Data::Vector{<:Real}; nbins::Int64 = 100, minValue::Real = Inf)
    m, n = extrema(Data)
    isinf(minValue) ? nothing : n > minValue ? m = minValue : error("minValue is greater than maximum value in data")
    
    bins = collect(range(m, stop = n, length = nbins))
    cDist = zeros(nbins)
    count = 1
    
    for value = bins
        cDist[count] = sum(Data .<= value)/length(Data)
        count = count + 1
    end
    
    return cDist, bins
end;

function get_dist(Data::Vector{<:Real}, bin_size::Float64; showNbins::Bool = false)
    min_value, max_value = extrema(Data)
    num_digits = convert(Int64, ceil(abs(log10(bin_size))))
    l = length(Data)
    
    number_of_bins = convert(Int64, ceil((max_value - min_value)/bin_size))
    distribution = zeros(number_of_bins + 1)
    
    for n = 1:number_of_bins
        bin_start = trunc(min_value + (n-1)*bin_size, digits = num_digits)
        bin_end = trunc(min_value + n*bin_size, digits = num_digits)
        
        distribution[n] = sum((Data .>= bin_start).*(Data .< bin_end))/l

        n == number_of_bins ? distribution[n+1] = sum(Data .== bin_end)/l : nothing
    end
    
    if showNbins
        return distribution, number_of_bins
    else
        return distribution
    end
end

function print_progress(mssg::String, current::Int64, termination::Int64)
    p = 100*current/termination
    p < 1.0 ? begin print("\r"); print(mssg, " |", " "^100, "| ") end : nothing

    try
        n = convert(Int64, p)
        print("\r")
        if n < 100 
            print(mssg, " |", "\u220e"^n, " "^(100 - n), "| ", string(n), "\uff05")
        else
            print(mssg, " |", "\u220e"^n, " "^(100 - n), "| ", string(n), "\uff05\n")
        end
    catch
        nothing
    end
end

model_sigmoid(x, p) = p[1]./(p[2] + p[3]*exp.(-p[4]*(x - p[5])))
model_poly(x, p) = p[1]*x.^3 + p[2]*x.^2 + p[3]*x + p[4]
model_rc(x, p) = p[1] - p[2]*exp.(p[3]*x)
