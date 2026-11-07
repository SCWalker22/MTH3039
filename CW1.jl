"""
Function for q1

Args:
    u: Initial conditions
    p: Params
    t: current_time

Returns:
    Vector(dx, dy, dz)
"""
function q1(u, p, t)
    x, y, z = u
    a, b, c = p
    dx = - y - z
    dy = x + a*y
    dz = b + z*(x-c)
    return SVector(dx, dy, dz)
end

"""
"""
function psos1(A, index, value)
    # Estimate poincare section of traj(A) with hyperplane defined for when var(traj) with given index
    # crosses a specific value (- => +)
    # Use interpolation
    prev_state = A[1]
    crossings = []
    for point in A
        current_state = point
        if prev_state[index] < value && current_state[index] > value
            push!(crossings, current_state)
        end
        prev_state = current_state
    end
    return crossings
end

"""
Function to calculate the crossing times of a trajectory over a hyperplane using the poincare map and interpolation

Args:
    A: ODE
    tvec: Vector of times (used as end-point of measurings)
    index: Index of variable to use, ie x = 1, y = 2
    value: Value of fixed variable on hyperplane

Returns:
    Vector(float): Crossing times up to (strictly less than) maximum(tvec)
"""
function psos2(A, tvec, index, value)
    prev_state = A[1]
    crossing_times = []
    for (i, point) in enumerate(A)
        current_state = point
        if prev_state[index] < value && current_state[index] > value
            current_time = tvec[i]
            push!(crossing_times, current_time)
        end
        prev_state = current_state
    end
    return crossing_times
end

"""
"""
function NP_DS(u, p, t)
    N, P = u
    r, δ = p
    c, μ, ν, α, β, χ = 0.19, 0.03, 0.003, 800, 1.5, 0.004

    s = α*N*P/(β + N)

    dN = r*N* (1 - (c/r)*N) * ((N - μ)/(N + ν)) - s
    dP = χ*s - δ*P

    return SVector(dN, dP)
end

function simple_grouping(features, r)
    groups = [features[1]]
    group_indeces = [1]
    for point in features[2:end]
        set_flag = false
        for (i, g) in enumerate(groups)
            # println(i)
            # println("g = $g")
            if norm(g - point) < r
                push!(group_indeces, i)
                set_flag = true
                break
            end
        end
        if !set_flag
            push!(groups, point) # Create a new point in the group
            push!(group_indeces, length(groups))#Set a new group
        end
    end
    return group_indeces
end

"""
"""
function fg1(ds, featuriser, r, ics; Ttr, Dt, T)
    features = []
    for u0 in ics
        X, t = trajectory(ds, T, collect(u0); Ttr, Dt)
        feature = featuriser(X)
        push!(features, feature)
    end
    groups = simple_grouping(features, r)
    return groups
end

"""
"""
function fg2(ds, featuriser, r, ics; Ttr, Dt, T)
    features = []
    for point in ics
        u0 = collect(point)
        X, t = trajectory(ds, T, u0; Ttr, Dt)
        feature = featuriser(X)
        push!(features, feature)
    end
    groups = simple_grouping(features, r)
    group_init_traj = []
    boa_fractions = []
    for i in 1:maximum(groups)
        index = findfirst(==(i), groups)
        init_cond = collect(ics[index])
        X, t = trajectory(ds, T, init_cond; Ttr, Dt)
        push!(group_init_traj, X)
        # Now calculate BoA fractions
        num = count(isequal(i), groups)
        frac = num/length(groups)
        push!(boa_fractions, frac)
    end

    return (groups, group_init_traj, boa_fractions)
end

"""
"""
function duffing_oscillator(u, p, t)
    x, v = u
    δ, f, ω = p
    dx = v
    dv = x - x^3 - δ*v + f*cos(ω*t)
    return SVector(dx, dv, t)
end

"""
"""
function pp_bmap(ds::DynamicalSystem, ics)
pow = 4
grid_recur = (
range(0, 15^(1/pow); length = 100).^pow,
range(0, (30e-3)^(1/pow); length = 100).^pow,
)
bmap_recur = BasinMapRecurrences(ds, grid_recur;
consecutive_recurrences = 4000, consecutive_attractor_steps = 10,
Δt = 0.2, stop_at_Δt = true, store_once_per_cell = false,
)
labels = map(bmap_recur, ics)
end