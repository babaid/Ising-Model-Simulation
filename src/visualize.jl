module Vis
export gif_from_data
using Plots

function gif_from_data(data)
    
    anim = @animate for plane in data
        heatmap(plane, axis=nothing, c=:grays)
    end
    anim_cntr = 0
    while isfile(string("../plots/", anim_cntr, "_anim_ising2d.gif"))
           anim_cntr+=1
    end
    
    gif(anim, string("../plots/", anim_cntr, "_anim_ising2d.gif"), fps=25) 
    
end
end
