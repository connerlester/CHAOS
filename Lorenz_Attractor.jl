
using DynamicalSystems,NearestNeighbors
using GLMakie
using LaTeXStrings

##### LORENZ SYSTEM PARAMETERS !!!!
Base.@kwdef mutable struct Lorenz 
    dt::Float64 = 3*10^-3 # time step 
    ρ::Float64  = 45. 
    σ::Float64  = 20.
    β::Float64  = 8/3
    T_trans::Float64  = 1. # transient time
    N_run::Int64 = 5*10^5 ## number o time steps
    T_run::Float64 = N_run*dt ## run time
    r₀_base::Vector = [1.,10.,0.] ## initial state
    r₀_near::Vector = [.95,10.3,.2] ## "..." perturbed
    n::Int64 = 1500 ## -- track neighbors for "n" time steps
    ds = Systems.lorenz(ρ=ρ,σ=σ,β=β) ## set system (DynamicalSystems.jl)
end


### TRAJECTORIES !!!!
function phase(system;p = system)  
    ## need two trajectories to analyze near neighbors with...
    base_state = trajectory(p.ds,p.T_run, p.r₀_base, Δt = p.dt, Ttr = p.T_trans); ## traj 1
    near_state = trajectory(p.ds,p.T_run, p.r₀_near, Δt = p.dt, Ttr = p.T_trans); ## traj near 1
    return base_state,near_state
end

## FIND BEST FRIENDS !!!
function neighbors( state , near_state )
    ## fast alg to find near neighbors (NearestNeighbors.jl)
    tree = KDTree(transpose(Matrix(near_state))) 
    gotcha = []
    for t ∈ 1:length(near_state)
        x,y,z = state[t,1],state[t,2],state[t,3] ## point in "base" state space
        buddy = nn(tree,[x,y,z])[1] ## nearest point in "near" state space
        push!(gotcha, buddy) ## gotcha, buddy! 
    end
    return gotcha 
end

### RECORD NEAREST NEIGHBORS DISPLACEMENT !!!
function Δneighbors( state , near_state , system ;p = system,n=p.n)

    friends = neighbors(state , near_state); ## get near neighbors
    
    T = length(state) ## total number o time steps

    Δfriends = zeros(n,T) # init
    
    for t ∈ 1:T-n

        t_nb = friends[t] ## nearest friend at time t

        if t_nb+n<T && t_nb-n>0 ## stay in time bounds...
            
            t⁺ , t_nb⁺= t:t+n-1 , t_nb:t_nb+n-1 ## time horizons (n)

            x,y,z = state[t⁺,1],state[t⁺,2],state[t⁺,3] ## "base" trajectory ∈ [t,t+n]
            i,j,k = near_state[t_nb⁺,1],near_state[t_nb⁺,2],near_state[t_nb⁺,3] ## "near" trajectory ∈ [t,t+n]

            dont_go = @. sqrt((x-i)^2 + (y-j)^2 + (z-k)^2); ## measure distance between phase space points
            Δfriends[:,t] = dont_go ## record

        end
    end
    return Δfriends
end

### RUN EM
function displacement(system)
    state , near_state = phase(system) ## get trajectories
    Δ = Δneighbors(state , near_state,system) ## track neighbor displacement
    return Δ,state
end


##### GET DAT DAT !!!
Δ,state = displacement(Lorenz());

x,y,z = state[:,1],state[:,2],state[:,3]; ## trajectory coords (to visualize attractor...)

####### ANIMATION !!!!
function CHAOS()

    ## init fig stuff
    set_theme!(backgroundcolor = :black) 
    fig = Figure(resolution=(1500,1100))
    ax = Axis3(fig[1,1],elevation=π/25,azimuth=π/1.325)
    hidedecorations!(ax)  # hides ticks, grid and lables
    hidespines!(ax) 

    t= Observable(1) ## init time-step observable to update
    Δₜ = @lift(Δ[$t,:]-Δ[1,:]) ## init relative-displacement-o-neighbors observable to update

    ## plot attractor with displacement as color
    butterfly=lines!(ax, x,y,z,  color=Δₜ#=<-to be updated=#,  colorrange = (0,15),linewidth=3 ,colormap = :gnuplot)

    # blah
    Colorbar(fig[2, 1][1,4],butterfly,label=L"Δ(t)-Δ(0)",
    labelsize=35,ticklabelsize=20,labelcolor="white",ticklabelcolor="white",tickcolor="white",
    vertical = false)

    ## record the action !!!
    record(fig, "Lorenz_Attractor.mp4", 1:Lorenz().n; ## iterate over time and update "t"
        framerate=60) do next_t
            t[] = next_t    
    end
    
end



#####
CHAOS() #  :D
#####