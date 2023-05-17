function secondsToTime(s)

    secs = floor(Int, s % 60)
    mins = floor(Int, s / 60.0)

    h = floor(Int,mins / 60)

    mins = mins % 60

    string("Time: ", h, ":", mins, ":", secs)

end

function simulationRound(grid, α, τ, h2)
    
    sze = size(grid)
    result = zeros(sze)

    # Dirlichlit
    result[1, :] .= 100.0

    for i in 2:size(grid,1)-1
        for j in 2:size(grid,2)-1
            result[i,j] = grid[i,j] + τ * α * (( 
                      grid[i + 1, j] 
                    + grid[i - 1, j] 
                    + grid[i, j + 1]
                    + grid[i, j - 1]
                    - 4 * grid[i, j])

                    / h2 )
        end
    end

    # Neumann
    result[sze[1], :] .= result[sze[1]-1, :]
    result[:, sze[2]] .= result[:, sze[2]-1]
    result[:, 1] .= result[:, 2]
end


function simulation()
    α = 3.8e-6 # m^2 / s
    temp_danger = 45.0 # °C
    length = 1.0 # m
    measure_danger = (400, 25) # idx

    grid_size = (500, 50) 

    flame_temp = 100.0 # °C
    flame_location = 1 # idx


    #τ = 0.01 # ?
    h = length / grid_size[1]
    h2 = h * h
    τ = h2/ (2 * α)

    ##

    stick = zeros(grid_size)

    stick[flame_location, :] .=  flame_temp;
    t = 0

    
    while stick[measure_danger[1], measure_danger[2]] < temp_danger
        stick .= simulationRound(stick, α, τ, h2)
        (t % 20000 == 0) && print(t,τ, stick[measure_danger[1], measure_danger[2]])
        t += 1
    end

    println(secondsToTime(t * τ))
end

function print(t, τ, temp) 
    println(secondsToTime(t * τ), " Temp: " , temp)
    false
end

simulation()
