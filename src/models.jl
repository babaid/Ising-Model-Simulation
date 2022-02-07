module Ising

export run, Lattice, Magnetization


using StatsBase, Plots, ProgressBars, IJulia


"""Inits a lattice, no boundary conditions implemented"""
function Lattice(size, cold=true)
    if cold
        s = sample([-1], size)
    else
        s = sample([1, -1], size)
    end
    
    return s
end


"""
Hamiltonian without external magnetic field
"""
function H_0(lattice, J)
    N = length(lattice[:, 1])
    M = length(lattice[1, :])
    E=0
    for i in 2:N-1
        for j in 2:M-1
            E += -J*lattice[i, j]*(lattice[i+1, j]+lattice[i-1, j]+lattice[i, j+1]+lattice[i, j-1])
        end
    end
    return E/2 #because of double counting
end


"""
Hamiltonian with external magnetic field
"""
function H_B(lattice, J, B)
    E=0
    N = length(lattice[:, 1])
    M = length(lattice[1, :])
    for i in 2:N-1
        for j in 2:M-1
            E += -J*lattice[i, j]*(lattice[i+1, j]+lattice[i-1, j]+lattice[i, j+1]+lattice[i, j-1])
        end
    end
    E*=0.5
    E += B* sum(lattice)
    return E 
end

"""Updates one lattice point"""
function step!(s, J, B, T)
    
    b = 1/(T*1.38e-23)
    
    #calculate energy
    e_1 = H_B(s, J, B)
    
    #flip a random spin
    i, j = rand(2:length(s[:, 1])), rand(2:length(s[1, :])) 
    s[i, j] = -s[i, j]
    
    #periodic boundaries
    s[:, 1] = s[:, end]
    s[1, :] = s[end, :]
    
    #calculate  new energy
    e_2 = H_B(s, J, B)
    
    #energy difference
    de = e_2-e_1
    
    #Metropolis Hastings
    
    if de<0
        return
    else   
        choice=sample([1, 0], ProbabilityWeights([exp(-b*de), 1-exp(-b*de)]))
    
        if choice==1
            return
        else
            #back to state before if rejected
            s[i, j] = -s[i, j]
            s[:, 1] = s[:, end]
            s[1, :] = s[end, :]
            return
        end
    end
end   
function step!(s, J, T)
    
    b = 1/(T*1.38e-23)
    
    #calculate energy
    e_1 = H_0(s, J)
    
    #flip a random spin
    i, j = rand(2:length(s[:, 1])), rand(2:length(s[1, :])) 
    s[i, j] = -s[i, j]
    
    #periodic boundaries
    s[:, 1] = s[:, end]
    s[1, :] = s[end, :]
    
    #calculate  new energy
    e_2 = H_0(s, J)
    
    #energy difference
    de = e_2-e_1
    
    #Metropolis Hastings
    
    if de<0
        return
    else   
        choice=sample([1, 0], ProbabilityWeights([exp(-b*de), 1-exp(-b*de)]))
    
        if choice==1
            return
        else
            #back to state before if rejected
            s[i, j] = -s[i, j]
            s[:, 1] = s[:, end]
            s[1, :] = s[end, :]
            return
        end
    end
end   




    
""" Runs simulation for a given configuration, saves an animation, returns magnetization"""
function run(config)
    T_start, T_end, B, J = 0, 0, 0, 0
    cold = true
    size = nothing
    #configuration
    for (key, value) in config
        if key=="T_start"
            T = value
        elseif key=="T_end" 
            T_end=value
        elseif key=="B" #magnetic field
            B = value
        elseif key=="J" #defines energy scale and critical temperature
            J = value
        elseif key == "cold_start"
            cold = value
        elseif key =="size"
            size = value
        
        end
        
    end
    
    #start the simulation
    loops = 3*300*300
    s = Lattice(size, cold)
    s[:, 1] = s[:, end]
    s[1, :] = s[end, :]
    #25 frames per second, 1000 frames in total, animation will be 40 seconds
    dT = (T_end -T_start)/1000
    iter = ProgressBar(T_start:dT:T_end)
    magnetization = []
    if B == 0
        anim = @animate for T in iter
            
           
            IJulia.clear_output(true)
            println(iter, "Temperature $T K")
            for j in 1:loops
                step!(s, J, T)
            end
            
            
            push!(magnetization, sum(lattice)/length(lattice))
                
    
            heatmap(s, axis=nothing, c=:grays, aspect_ratio=1)
        end
        
    else
        anim = @animate for T in iter
            
            
            IJulia.clear_output(true)
            println(iter, "Temperature $T K")
            for j in 1:loops
                step!(s, J, B, T)
            end
            
            push!(magnetization, sum(lattice)/length(lattice))
                
            
            heatmap(s, axis=nothing, c=:grays, aspect_ratio=1)
        end
    end   
    
    gif(anim, "../plots/simulation2.gif", fps=25)
    
    return magnetization
    
end

function Magnetization(s)
    return sum(s)/length(s)
end



end