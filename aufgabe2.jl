function secondsToTime(s)

    secs = floor(Int, s % 60)
    mins = floor(Int, s / 60.0)

    h = floor(Int,mins / 60)

    mins = mins % 60

    string("Time: ", h, ":", mins, ":", secs)

end

function simulationRound(grid, α, τ, h)
    
    sze = size(grid)
    result = zeros(sze)

    result[1] = 100.0

    for i in 2:size(grid,1)-1
        result[i] = grid[i] + τ * α * (( grid[i + 1] + grid[i - 1] - 2 * grid[i] ) /  (h * h))
    end
    result[sze[1]] = result[sze[1]-1]
    result
end


function simulation()
    α = 3.8e-6 # m^2 / s
    temp_danger = 45.0 # °C
    length = 1.0 # m
    measure_danger = 800 # idx

    grid_size = (1000) 

    flame_temp = 100.0 # °C
    flame_location = 1 # idx


    τ = 0.01 # ?
    h = length / grid_size[1]

    ##

    stick = zeros(grid_size)

    stick[flame_location] =  flame_temp;
    t = 0
    while stick[measure_danger] < temp_danger
        stick .= simulationRound(stick, α, τ, h)
        t += 1
    end

    println(secondsToTime(t * τ))
end


simulation()
