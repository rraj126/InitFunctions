using PyPlot, PyCall
using StatsBase, Statistics, LinearAlgebra

export bandplot, gridplot, gridplot2, histplot, bxplot, vlnplot, cdistplot, heatplot, barplot, fitplot

ArrayOfArray = Union{Array{Array{Float64, 1}, 1}, Array{Array{Int64, 1}, 1}}

function get_plot_index(row_config::Vector{Int64})
    plot_index = zeros(Int64, sum(row_config))
    j = 1
    for i in 1:length(row_config)
        init_index = (i-1)*maximum(row_config)
        if !iszero(row_config[i])
            plot_index[j:j+row_config[i]-1] = init_index .+ collect(1:row_config[i])
            j = j + row_config[i] 
        
        end
    
    end

    return plot_index
end

function bandplot(Data::Union{ArrayOfArray, Array{<:Real}}; showMeans::Bool = true, band::String = "", color = "", mean_color = "", match_color::Bool = false, label::String = "")
    data_type, number_of_boxes = get_nboxes(Data)
    
    data_type == "array" && number_of_boxes == 1 ? error("cannot plot single vector") : nothing
    y_error = zeros(number_of_boxes, 2)
    y_values = zeros(number_of_boxes)
    data_means = zeros(number_of_boxes)
    
    isempty(band) ? band = "std" : nothing
    band == "std" || band == "iqr" ? nothing : error("choose band to represent mean and standard deviation (std) or median and inter-quartile range (iqr)")
    
    for boxNo = 1:number_of_boxes
        data_type == "array" ? data = Data[boxNo] : data_type == "matrix" ? begin data = @view Data[:, boxNo] end : error("data type not determined")
        data_means[boxNo] = mean(data)
        
        if band == "iqr"
            temp = summarystats(data)
            y_values[boxNo] = median(data)
            y_error[boxNo, 1] = temp.q25
            y_error[boxNo, 2] = temp.q75
        else
            y_values[boxNo] = mean(data)
            y_error[boxNo, 1] = mean(data) + std(data)
            y_error[boxNo, 2] = mean(data) - std(data)
        end
    end
    
    ax = gca()
    ax.set_title("dp", fontsize = 0)
    
    isempty(color) ? color = COLOR_CYCLE[1] : nothing
    ax.fill_between(collect(1:number_of_boxes), vec(y_error[:, 1]), vec(y_error[:, 2]), alpha = 0.5, color = color, linewidth = 0.0)
    
    if isempty(label) 
        ax.plot(collect(1:number_of_boxes), y_values, color = color, linewidth = 3)
    else
        ax.plot(collect(1:number_of_boxes), y_values, color = color, linewidth = 3, label = label)
    end
    
    match_color ? mean_color = color : isempty(mean_color) ? mean_color = COLOR_CYCLE[2] : nothing
    showMeans && band == "iqr" ? ax.scatter(collect(1:number_of_boxes), data_means, color = mean_color, s = 150, marker = "d", zorder = 3, alpha = 0.5) : nothing
end;
                                                    
function gridplot(X::Array{<:Real, 2}; atom_dims::Tuple{Int64, Int64} = (0, 0), cmap = "magma", trans::Bool = false, sample::Bool = false, reverse::Bool = true, scale_by_matrix::Bool = false, scale_by_column::Bool = false)
    dim1, dim2 = size(X)
    sample ? number_of_atoms = min(50, dim2) : number_of_atoms = dim2
    temp = auto_reshape(Array{Float64, 1}(undef, number_of_atoms))
    grid_rows, grid_cols = size(temp)
    
    a_max = zero(eltype(X))
    a_min = zero(eltype(X))
    scale_by_matrix && scale_by_column ? scale_by_matrix = false : nothing

    fig = figure(figsize = (grid_cols, grid_rows))
    fig.patch.set_facecolor("#000000")
    
    for atomNo = 1:number_of_atoms
        sum(atom_dims) == 0 ? atom = auto_reshape(vec(X[:, atomNo])) : atom = reshape(vec(X[:, atomNo]), atom_dims[1], atom_dims[2])
        trans ? atom = atom' : nothing
        reverse ? atom = reverse_matrix(atom) : nothing
        
        scale_by_matrix ? a = maximum(abs.(X)) : scale_by_column ? a = maximum(abs.(atom)) : begin a = Inf; a_max = maximum(atom); a_min = minimum(atom) end
        isinf(a) ? nothing : begin a_max = a; a_min = -a end 
        
        ax = fig.add_subplot(grid_rows, grid_cols, atomNo)
        pcolormesh(atom, cmap = cmap, vmin = a_min, vmax = a_max)
        xticks([])
        yticks([])
        set_axprops(ax, false)
    end
end;

function histplot(Data::Array{<:Real, 1}; bins::Array{<:Real, 1} = zeros(1), normalized::Bool = true, color = "", alpha::Float64 = 1.0)
    ax = gca()
    isempty(color) ? color = COLOR_CYCLE[1] : nothing                                                    
    sum(bins) == 0 ? (values, bins) = ax.hist(Data, color = color, alpha = alpha) : (values, bins) = ax.hist(Data, bins, color = color, alpha = alpha)
    
    maxValue_normalized = maximum(values)/length(Data)
    num_digits = ceil(abs(log10(maxValue_normalized)))
    maxValue_normalized_ceil = ceil(maxValue_normalized, digits = convert(Int64, num_digits))

    y_pos = collect(0.0:10^(-num_digits):maxValue_normalized_ceil)
    y_values = y_pos*length(Data)

    normalized ? yticks(y_values, y_pos) : nothing
end;

function bxplot(Data::Union{ArrayOfArray, Array{<:Real}}; pos::Array{Int64, 1} = zeros(Int64, 1), showMeans::Bool = true, showFliers::Bool = true, color = "", label::String = "", swarm::Bool = false)
    ax = gca()
    data_type, number_of_boxes = get_nboxes(Data)

    ws = 0.7*ones(number_of_boxes)
    sum(pos) == 0 ? pos = collect(1:number_of_boxes) : nothing 
    bp = ax.boxplot(Data, widths = ws, positions = pos, showmeans = showMeans, showfliers = showFliers, patch_artist = true)
    
    ax.set_title("dp", fontsize = 0)
    isempty(color) ? color = COLOR_CYCLE[1] : nothing
    isempty(label) ? config_bxplot(bp, color) : config_bxplot(bp, color, label = label)
    
    if swarm 
        for boxNo = 1:number_of_boxes
            data_type == "matrix" ? box_data = vec(Data[:, boxNo]) : nothing
            if data_type == "array"
                number_of_boxes == 1 ? box_data = vec(Data) : number_of_boxes > 1 ? box_data = vec(Data[boxNo]) : error("No data to swarm")
            end
            
            xvalues = vec(pos[boxNo] .+ (0.5*rand(length(box_data)) .- 0.25))
            ax.scatter(xvalues, box_data, marker = ".", alpha = 0.5, color = color, s = 80, zorder = 1)
        end
    end      
end;

function vlnplot(Data::Union{ArrayOfArray, Array{<:Real}}; pos::Array{Int64, 1} = zeros(Int64, 1), showMeans::Bool = true, showMedians::Bool = true, color = "", label::String = "", swarm::Bool = false)
    ax = gca()
    data_type, number_of_boxes = get_nboxes(Data)

    ws = 0.7*ones(number_of_boxes)
    sum(pos) == 0 ? pos = collect(1:number_of_boxes) : nothing 
    vp = ax.violinplot(Data, widths = ws, positions = pos, showmeans = showMeans, showmedians = showMedians)
    
    ax.set_title("dp", fontsize = 0)
    isempty(color) ? color = COLOR_CYCLE[1] : nothing
    isempty(label) ? config_vlnplot(vp, color) : config_vlnplot(vp, color, label = label)
    
    if swarm 
        for boxNo = 1:number_of_boxes
            data_type == "matrix" ? box_data = vec(Data[:, boxNo]) : nothing
            if data_type == "array"
                number_of_boxes == 1 ? box_data = vec(Data) : number_of_boxes > 1 ? box_data = vec(Data[boxNo]) : error("No data to swarm")
            end
            
            xvalues = vec(pos[boxNo] .+ (0.5*rand(length(box_data)) .- 0.25))
            ax.scatter(xvalues, box_data, marker = ".", alpha = 0.5, color = color, s = 80, zorder = 1)
        end
    end
end;

function cdistplot(Data::Array{<:Real, 1}; nbins::Int64 = 100, style::String = "line", color = "", minValue::Real = Inf)
    cDist, bins = get_cdist(Data, nbins = nbins, minValue = minValue)
    
    ax = gca()
    ax.set_ylim(0, 1)
                                                        
    if isempty(color)
        style == "bar" ? ax.bar(bins, cDist, align = "center") : style == "line" ? ax.plot(bins, cDist, "-") : nothing
    else
        style == "bar" ? ax.bar(bins, cDist, align = "center", color = color) : style == "line" ? ax.plot(bins, cDist, "-", color = color) : nothing
    end
end;

function heatplot(Data::Array{<:Real, 2}; contours::Bool = false, levels::Array{<:Real, 1} = Array{Int64, 1}(undef, 0), cmap = "", contour_color = "", vmin = 0, vmax = 0)
    ax = gca()
    
    m, n = size(Data)
    y_mesh = repeat(collect(0:m), 1, n + 1)
    x_mesh = repeat(collect(0:n)', m + 1, 1)

    #isa(cmap, ColorMap) ? nothing : cmap = "magma"
    isempty(cmap) ? cmap = "magma" : nothing
    
    vmin + vmax == 0 ? pcolormesh(x_mesh, y_mesh, Data, cmap = cmap) : pcolormesh(x_mesh, y_mesh, Data, cmap = cmap, vmin = vmin, vmax = vmax)
    
    if contours
        y_contours = repeat(collect(1:m), 1, n) .- 0.5
        x_contours = repeat(collect(1:n)', m, 1) .- 0.5
        
        isempty(contour_color) ? contour_color = COLOR_CYCLE[1] : nothing
        isempty(levels) ? begin data_min, data_max = extrema(vec(Data)); levels = collect(range(data_min, stop = data_max, length = 10)) end : nothing
        int_levels = levels[rem.(levels, 1) .== 0]
        flt_levels = levels[rem.(levels, 1) .!= 0]

        cs = ax.contour(x_contours, y_contours, Data, levels = levels, colors = contour_color)
        ax.clabel(cs, int_levels, fmt = "%0.0f")
        ax.clabel(cs, flt_levels, fmt = "%0.2f")
    end
end;

heatplot(Data::BitArray{2}) = heatplot(1 .* Data)
                                        
function barplot(Data::Union{ArrayOfArray, Array{<:Real}}; pos::Array{<:Real, 1} = zeros(Float64, 1), yerror::Bool = true, bottom = "", color = "", label::String = "")
    ax = gca()
    
    data_type, number_of_boxes = get_nboxes(Data)
    isempty(color) ? color = COLOR_CYCLE[1] : nothing
    
    if number_of_boxes == 1
        sum(pos) == 0 ? pos = collect(1:length(Data)) : nothing
        isempty(label) ? ax.bar(pos, Data, align = "center", width = 0.7, color = color) : ax.bar(pos, Data, align = "center", width = 0.7, color = color, label = label)
    else
        if data_type == "array"
            sum(pos) == 0 ? pos = collect(1:length(Data)) : nothing
            
            mean_data = zeros(length(Data))
            std_data = zeros(length(Data))
            for arrayNo = 1:length(Data) begin mean_data[arrayNo] = mean(Data[arrayNo]); std_data[arrayNo] = std(Data[arrayNo]) end end
                                                        
        elseif data_type == "matrix"
            sum(pos) == 0 ? pos = collect(1:size(Data, 2)) : nothing
            
            mean_data = vec(mean(Data, dims = 1))
            std_data = vec(std(Data, dims = 1))
        
        else
            error("Data type cannot be determined")
        
        end
        
        l = maximum(pos) - minimum(pos) + 1
        if yerror
            if isempty(label)
                ax.bar(pos, mean_data, align = "center", width = 0.7, color = color, yerr = std_data, capsize = 3000/(32*l))
            else
                ax.bar(pos, mean_data, align = "center", width = 0.7, color = color, yerr = std_data, capsize = 3000/(32*l), label = label)
            end
        else
            if isempty(label)
                ax.bar(pos, mean_data, align = "center", width = 0.7, color = color)
            else
                ax.bar(pos, mean_data, align = "center", width = 0.7, color = color, label = label)
            end
        end
    end
end

function fitplot(T::DataType, x::Vector{<:Real}, y::Vector{<:Real}; label::String = "", ptype::String = "scatter")
    ax = gca()
    if ptype == "scatter"
        isempty(label) ? ax.scatter(x, y) : ax.scatter(x, y, label = label)
    end

    xfit, yfit = fit_curve(T, x, y)
    isempty(label) ? ax.plot(xfit, yfit) : ax.plot(xfit, yfit, label = "fit")

end


function gridplot2(data::Union{Matrix{<:Real}, Dict{Int64, Matrix{<:Real}}}, row_config::Vector{Int64};  atom_dims::Tuple{Int64, Int64} = (0, 0), cmap = "magma", trans::Bool = false, reverse::Bool = true)
    grid_rows, grid_cols = length(row_config), maximum(row_config)

    fig = figure(figsize = (grid_cols, grid_rows))
    fig.patch.set_facecolor("#000000")

    plot_index = get_plot_index(row_config)
    
    isa(data, Matrix) ? number_of_atoms = size(data, 2) : number_of_atoms = length(data)
    number_of_atoms != length(plot_index) ? error("data cannot be plotted in specified format") : nothing
    
    for atomNo in 1:number_of_atoms
        if isa(data, Matrix)
            if iszero(sum(atom_dims)) 
                atom = auto_reshape(data[:, atomNo])

            else
                atom = reshape(vec(X[:, atomNo]), atom_dims[1], atom_dims[2])
            end
            
        else
            atom = data[atomNo]
        end

        trans ? atom = atom' : nothing
        reverse ? atom = reverse_matrix(atom) : nothing
        
        ax = fig.add_subplot(grid_rows, grid_cols, plot_index[atomNo])
        pcolormesh(atom, cmap = cmap)
        xticks([])
        yticks([])
        set_axprops(ax, false)
    end
end
