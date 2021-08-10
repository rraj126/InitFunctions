using PyCall, PyPlot

export my_cmap

sns = pyimport("seaborn")
#clrs = pyimport("matplotlib.colors")

c = sns.color_palette("rocket", 256)
c_array = Array{Float64, 2}(undef, 256, 3)
for i in 1:256
    c_array[i, 1], c_array[i, 2], c_array[i, 3] = c[i]
end

my_cmap = ColorMap("rocket", c_array)