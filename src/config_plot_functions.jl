using PyPlot
using PyCall

COLOR_CYCLE = ["#00bfff"; "#ff0000"; "#69ef69"; "#cc3fff"; "#ff7c6d"]
export create_axes, plotNo, config_plots, COLOR_CYCLE


function create_axes(N::Int64; projections::Union{String, Vector{String}} = "")
    projection_array = Vector{String}(undef, N)
    if isa(projections, String)
        isempty(projections) ? projection_array .= "rectilinear" : projection_array .= projections
    else
        projection_array = repeat(projections, ceil(Int64, N/length(projections))) 
    end

    m, n = size(auto_reshape(trues(N)))
    axes = Matrix{PyObject}(undef, m, n); ax_count = 1

    for i in 1:m, j in 1:n
        axes[i, j] = subplot(m, n, ax_count, projection = projection_array[ax_count])
        ax_count += 1
    end

    fig = gcf(); fig.set_size_inches(6.5*n, 4*m)
    
    return axes
end;


function check_heatmap(ax)
    ax._current_image == nothing ? presence = false : ax._current_image.__class__.__name__ == "QuadMesh" ? presence = true : presence = false
    
    return presence
end;


function check_distplot(ax)
    t = ax.get_title()
    t == "dp" ? bp_present = true : bp_present = false
    
    return bp_present
end;


function plot_colorbar(fig, ax)
    ax_pos = ax.get_position()
    
    x1 = ax_pos.x1
    x0 = ax_pos.x0
    y0 = ax_pos.y0
    w = ax_pos.width
    h = ax_pos.height
    
    cax_pos = (x0 + 0.97*w, y0, 0.03*w, h)
    ax_pos_new = (x0, y0, 0.95*w, h)
    
    ax.set_position(ax_pos_new)
    cax = fig.add_axes(cax_pos)
    cb = fig.colorbar(ax._current_image, cax = cax)
    
    return cb
end;


function set_axprops(ax, legends::Bool; c = "darkgray", cb_handle = nothing, frameon::Bool = false)   
    ax.set_facecolor("#000000")
    ax.set_axisbelow(true)

    x_label = ax.get_xlabel()
    y_label = ax.get_ylabel()
    
    polar = haskey(ax.spines, "polar")
    
    if frameon
        if polar
            ax.spines["polar"].set_color(c)
        else
            ax.spines["bottom"].set_color(c)
            ax.spines["top"].set_color(c)
            ax.spines["left"].set_color(c)
            ax.spines["right"].set_color(c)    
        end
    else
        ax.grid(linestyle = "--", color = "darkgray", alpha = 0.3)
        if polar
            ax.spines["polar"].set_visible(false)
        else
            ax.spines["top"].set_visible(false)
            ax.spines["bottom"].set_visible(false)
            ax.spines["left"].set_visible(false)
            ax.spines["right"].set_visible(false)
        end
    end

    polar ? ax.tick_params(axis = "x", colors = c, labelsize = 8) : ax.tick_params(axis = "x", colors = c, labelsize = 8, rotation = 45)
    ax.tick_params(axis = "y", colors = c, labelsize = 8)

    for tick in ax.get_xticklabels() begin tick.set_fontname("Courier New"); tick.set_fontweight("bold") end end
    for tick in ax.get_yticklabels() begin tick.set_fontname("Courier New"); tick.set_fontweight("bold") end end

    ax.set_xlabel(x_label, fontsize = 10, fontname = "Courier New", fontweight = "bold", color = c)
    ax.set_ylabel(y_label, fontsize = 10, fontname = "Courier New", fontweight = "bold", color = c)
    
    if cb_handle != nothing
        ax.grid("off")
        cb_handle.outline.set_visible(false)
        cb_handle.ax.tick_params(labelsize = 8, colors = c)
        for tick in cb_handle.ax.get_yticklabels() begin tick.set_fontname("Courier New"); tick.set_fontweight("bold") end end
    end
    legends ? begin l = ax.legend(loc = 2, bbox_to_anchor = (1.02, 1.0), frameon = false); setp(l.get_texts(), color = c, fontname = "Courier New", fontweight = "bold", fontsize = 10); end : nothing
end;


function set_pltcolors(ax)
    line_list = ax.get_lines()
    patch_list = ax.patches
    collection_list = ax.collections
        
    !isempty(line_list) && !isempty(patch_list) && !isempty(collection_list) ? bar_with_error = true : bar_with_error = false
    
    if !bar_with_error
        if !isempty(line_list)
            for lineNo = 1:length(line_list) 
                index = convert(Int64, rem(lineNo, length(COLOR_CYCLE)))
                index == 0 ? clr_index = length(COLOR_CYCLE) : clr_index = index

                line_list[lineNo].set_color(COLOR_CYCLE[clr_index])
                line_list[lineNo].set_linewidth(3)
            end
        end

        if !isempty(collection_list)
            for collectionNo = 1:length(collection_list) 
                index = convert(Int64, rem(collectionNo, length(COLOR_CYCLE)))
                index == 0 ? clr_index = length(COLOR_CYCLE) : clr_index = index

                collection_list[collectionNo].set_color(COLOR_CYCLE[clr_index])
                collection_list[collectionNo].set_alpha(0.5)
                collection_list[collectionNo].set_sizes((150, ))
            end
        end
        !isempty(patch_list) ? for patchNo = 1:length(patch_list) patch_list[patchNo].set_color(COLOR_CYCLE[1]) end : nothing
    else
        for lineNo = 1:length(line_list) begin line_list[lineNo].set_color("gray"); line_list[lineNo].set_markeredgewidth(2) end end
        for collectionNo = 1:length(collection_list) begin collection_list[collectionNo].set_edgecolors("gray"); collection_list[collectionNo].set_linewidth(2) end end
    end
end;


function config_plots(axs...)
    fig = gcf()
    fig.patch.set_facecolor("#000000")
    
    plotted_axes = fig.get_axes()
    length(plotted_axes) == 1 ? fig.set_size_inches(6.5, 4) : nothing
    
    all_axes = collect(1:length(plotted_axes))
    
    if isempty(collect(axs))
        color_axes = all_axes
    else
        temp = collect(axs)
        temp isa Array{Int64, 1} ? inp_axes = temp : temp[1] isa Array ? inp_axes = temp[1] : temp[1] isa Tuple ? inp_axes = collect(temp[1]) : nothing

        inp_axes = vec(inp_axes)
        cpos = inp_axes[inp_axes .> 0]
        cneg = abs.(inp_axes[inp_axes .< 0])
        
        
        if !isempty(cpos) && isempty(cneg)
            color_axes = cpos
        
        elseif isempty(cpos) && !isempty(cneg)
            color_axes = setdiff(all_axes, cneg)
        
        else 
            error("Cannot accept both positive and negative axes simultaneously")
        end        
    end        
    
    for axisNo = 1:length(plotted_axes)
        ax = plotted_axes[axisNo]
        
        h, l = ax.get_legend_handles_labels()
        lp = !isempty(l)
        
        in_color_axes = axisNo in color_axes
        dp = check_distplot(ax)
        
        if check_heatmap(ax) 
            cb = plot_colorbar(fig, ax)
            set_axprops(ax, lp, cb_handle = cb) 
        
        elseif  in_color_axes && !dp 
            set_pltcolors(ax) 
            set_axprops(ax, lp)
        
        else
            set_axprops(ax, lp)
        
        end
    end
end;


function config_bxplot(bp, c; label = "")
    number_of_boxes = length(bp["boxes"])
    
    length(bp["means"]) > 0 ? means_present = true : means_present = false
    length(bp["fliers"]) > 0 ? fliers_present = true : fliers_present = false
    
    for index = 1:number_of_boxes
  
        setp(bp["boxes"][index], color = c, facecolor = c, alpha = 0.5)
        setp(bp["whiskers"][2*index - 1], color = c, linestyle = "-")
        setp(bp["whiskers"][2*index], color = c, linestyle = "-")
        setp(bp["caps"][2*index - 1], color = c)
        setp(bp["caps"][2*index], color = c)
        setp(bp["medians"][index], color = c)
        
        fliers_present ? setp(bp["fliers"][index], marker = "+", markeredgecolor = "gray") : nothing
        means_present ? setp(bp["means"][index], marker = "d", markeredgecolor = c, markerfacecolor = c, markersize = 8) : nothing
    end
    ax = gca()
                                        
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    
    ax.set_xlim(0, ceil(x2))
    isempty(label) ? nothing : ax.plot([0.5*x1 + 0.5*x2], [0.5*y1 + 0.5*y2], color = c, linewidth = 10, alpha = 0.5, label = label)
end;


function config_vlnplot(vp, c; label = "")
    if haskey(vp, "cmeans") 
        number_of_plots = length(vp["cmeans"].get_paths())
        mean_x = zeros(number_of_plots)
        mean_y = zeros(number_of_plots)
        count = 1
        
        for p in vp["cmeans"].get_paths()
            coord = p.vertices
            mean_x[count] = mean(coord[:, 1])
            mean_y[count] = mean(coord[:, 2])
            count = count + 1
        end
        
        ax = gca()
        ax.scatter(mean_x, mean_y, marker = "d", color = c, s = 100)
    end
    
    for key in keys(vp)
        part = vp[key]
        key == "cmeans" ?  part.set_visible(false) : key == "bodies" ? nothing : part.set_edgecolor(c)
    end
    
    for body in vp["bodies"]
        body.set_color(c)
        body.set_alpha(0.5)
        body.set_edgecolor("none")
    end
    
    ax = gca()

    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    isempty(label) ? nothing : ax.plot([0.5*x1 + 0.5*x2], [0.5*y1 + 0.5*y2], color = c, linewidth = 10, alpha = 0.5, label = label)
    
    ax.set_xlim(0, ceil(x2) + 0.3)
end;

function get_nboxes(Data::Union{Array{Array{Int64, 1}, 1}, Array{Array{Float64, 1}, 1}, Array{<:Real}})

    if isa(Data, Array{<:Real, 2})
        data_type = "matrix"
        size(Data, 1) != 1 ? number_of_boxes = size(Data, 2) : error("cannot plot a row vector")

    elseif isa(Data, Array{<:Real, 1})
        data_type = "array"
        number_of_boxes = 1

    else
        data_type = "array"
        1 in length.(Data) ? error("cannot plot data") : number_of_boxes = length(Data)

    end
    
    return data_type, number_of_boxes
end;
    

function plotNo(N...)
    
    length(N) == 0 ? axisNo = 1 : (axisNo,) = N
    
    fig = gcf()
    plotted_axes = fig.get_axes()
    number_of_axes = length(plotted_axes)

    m, n = size(auto_reshape(ones(number_of_axes)))
    number_of_axes < axisNo ? error("Axis doesn't exist") : subplot(m, n, axisNo)
end;


function polarplot()
    fig = gcf()
    plotted_axes = fig.get_axes()

    ax = gca()
    axis_position = ax.get_position()
    axis_index = findall(plotted_axes .== ax)
    ax.set_xticks([])
    ax.set_yticks([])

    polar_ax = fig.add_axes(axis_position, polar = true)
    polar_ax.set_yticks([])
end
