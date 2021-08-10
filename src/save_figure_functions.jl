using PyPlot
using PyCall

export save_figure

function get_figprops(fig)
    
    plotted_axes = fig.get_axes()
    number_of_axes = length(plotted_axes)
    heatmap_axes = falses(number_of_axes)
    legend_axes = falses(number_of_axes)
    
    for axisNo = 1:number_of_axes
        heatmap_axes[axisNo] = check_heatmap(plotted_axes[axisNo])
        h, l = plotted_axes[axisNo].get_legend_handles_labels()
        legend_axes[axisNo] = !isempty(l)
    end
        
    
    true_axes = number_of_axes - sum(convert(Array{Int64, 1}, heatmap_axes))
    cbars = heatmap_axes[1:true_axes]
    legends = legend_axes[1:true_axes]
    
    return cbars, legends, true_axes
end;

function save_figure(fname; fmt = "png", frameon::BitArray{1} = falses(1), config_axes = "", plot_type = "", invert_background = false)
    if(!isdir("Figures"))
        mkdir("Figures/")
    end
    
    filename = string("Figures/", fname, ".", fmt)
    
    fig = gcf()
    if invert_background
        if isempty(plot_type)
            cbars, legends, true_axes = get_figprops(fig)
            plotted_axes = fig.get_axes()
            figure_axes = plotted_axes[1:true_axes]

            length(frameon) == length(figure_axes) ? nothing : frameon = falses(length(figure_axes))

            for axisNo = length(figure_axes)+1:length(plotted_axes)
                fig.delaxes(plotted_axes[axisNo])
            end

            for axisNo = 1:length(figure_axes)
                ax = figure_axes[axisNo]
                lp = legends[axisNo]
                cp = cbars[axisNo]

                if cp
                    cb = plot_colorbar(fig, ax)
                    set_axprops(ax, lp, cb_handle = cb, c = "black", frameon = frameon[axisNo])
                else
                    set_axprops(ax, lp, c = "black", frameon = frameon[axisNo])
                end
            end

            savefig(filename, facecolor = "w", transparent = true, format = fmt, bbox_inches = "tight")

            new_axes = fig.get_axes()
            for axisNo = length(figure_axes)+1:length(new_axes) 
                fig.delaxes(new_axes[axisNo])
            end

            isempty(config_axes) ? config_plots() : begin axes = collect(config_axes); config_plots(axes) end
        end
        if plot_type == "gridplot"
            grid_axes = fig.get_axes()
            length(frameon) == length(grid_axes) ? nothing : frameon = falses(length(grid_axes))

            for axisNo = 1:length(grid_axes) begin ax = grid_axes[axisNo]; set_axprops(ax, false, c = "black", frameon = frameon[axisNo]) end end 
            savefig(filename, facecolor = "w", transparent = true, format = fmt)
            for axisNo = 1:length(grid_axes) begin ax = grid_axes[axisNo]; set_axprops(ax, false, frameon = frameon[axisNo]) end end 
        end   
    else
        savefig(filename, format = fmt, bbox_inches = "tight")
    end
end;
