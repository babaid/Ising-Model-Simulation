module Ising

export run, Lattice, Magnetization
include("visualize.jl")


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

"""Updates one lattice point, with periodic boundary conditions"""
function step_periodic!(s, J, B, T)
    
    b = 1/(T*1.38e-23)
    
    #calculate energy
    e_1 = H_B(s, J, B)
    
    #flip a random spin
    i, j = rand(2:length(s[:, 1])-1), rand(2:length(s[1, :])-1)
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


"""
One update step without periodic boundary conditions
"""
function step_non_periodic!(s, J, B, T)
    
    b = 1/(T*1.38e-23)
    
    #calculate energy
    
    e_1 = H_B(s, J, B)
    
    #flip a random spin
    i, j = rand(1:length(s[:, 1])), rand(1:length(s[1,:]))
    s[i, j] = -s[i, j]
    
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
            return
        end
    end
end   




""" Runs simulation for a given configuration, saves an animation, returns magnetization"""
function run(config::Dict{String, Any}=Dict([("T_start", 0), ("T_end", 500),("dT", 1), ("J", 2e-21), ("B", 0), ("size", (50, 50)), ("samples", 1000), ("periodic", false), ("cold_start", true)]) )
    
    conf = Dict([("T_start", 0), ("T_end", 500),("dT", 1), ("J", 2e-21), ("B", 0), ("size", (50, 50)), ("samples", 1000), ("periodic", false), ("cold_start", true)])
    
    
    
    #update configuration values
    for (key, value) in config
       conf[key]=value
    end
    
    
    
    #start the simulation
    #25 frames per second, 1000 frames in total, animation will be 40 seconds
    iter = ProgressBar(conf["T_start"]:conf["dT"]:conf["T_end"])
    magnetization = []
    spins =[]
    #new simulation logfile
    sim_cntr = 1
    
    
    s = Lattice(conf["size"], conf["cold_start"])
    
    if conf["periodic"]
        
        s[:, 1] = s[:, end]
        s[1, :] = s[end, :]
        
        for T in iter
            IJulia.clear_output(true)
            println(iter, "Temperature $T K")
            #the situation at a specific T
            
            for j in 1:conf["samples"]
                step_periodic!(s, conf["B"], conf["J"], T)
            end
            
             #save spins and magnetization
            push!(magnetization, sum(s)/length(s))
            push!(spins, copy(s))  
        end
       
        
    elseif !conf["periodic"]
        
        for T in iter
            
            IJulia.clear_output(true)
            println(iter, "Temperature $T K")
            
            for j in 1:conf["samples"]
                step_non_periodic!(s, conf["B"], conf["J"], T)
            end
            
            push!(magnetization, sum(s)/length(s))
            push!(spins, copy(s))  
        
            
        end

        
    end
    
    
    
    
    
    
    #Saving log files
    #save logfile
    while isfile(string("../data/sims/", sim_cntr, "_sim_data_ising2d"))
        sim_cntr+=1
    end
    
    spin_log = open(string("../data/sims/", sim_cntr, "_sim_data_ising2d.txt"), "a")
    #mag_log = open(string("../data/sims/", sim_cntr, "_sim_data_ising2d.txt"), "a")
    
    for plane in spins
        for  elem in plane
            
            write(spin_log, string(elem, "\t"))
              
        end
        write(spin_log, "\n")
    end
    
    close(spin_log)
    #close(mag_log)
    
    return (spins, magnetization)
    
end

function Magnetization(s)
    return sum(s)/length(s)
end




end