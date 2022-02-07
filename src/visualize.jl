module Ising.Vis

export anim

using Plots

function show_state(s)
    heatmap(s, axis=nothing, c:grays)
end



end
