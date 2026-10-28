function q1(u, p, t)
    x, y, z = u
    a, b, c = p
    dx = - y - z
    dy = x + a*y
    dz = b + z*(x-c)
    return SVector(dx, dy, dz)
end

function psos1(A, index, value)
    # Estimate poincare section of traj(A) with hyperplane defined for when var(traj) with given index
    # crosses a specific value (- => +)
    # Use interpolation
    plane = (index, value)
    pmap = poincaresos(A, plane, 1000, direction = -1)
    return pmap
end

"""
Function to calculate the crossing times of a trajectory over a hyperplane using the poincare poincare_map

Args:
    A: ODE
    tvec: Vector of times (used as end-point of measurings)
    index: Index of variable to use, ie x = 1, y = 2
    value: Value of fixed variable on hyperplane

Returns:
    Vector(float): Crossing times up to (strictly less than) maximum(tvec)
"""
function psos2(A, tvec, index, value)
    max_time = maximum(tvec)
    times = []
    plane = (index, value)
    pmap = poincaremap(A, plane, max_time, direction = -1)
    while true
        step!(pmap)
        t_cross = current_crossing_time(pmap)
        if t_cross < max_time
            push!(times, t_cross)
        else
            return times
        end
    end
end

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
        println("hi $i")
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

function duffing_oscillator(u, p, t)
    x, v = u
    δ, f, ω = p
    dx = v
    dv = x - x^3 - δ*v + f*cos(ω*t)
    return SVector(dx, dv, t)
end