using LinearAlgebra

function secondsToTime(s)

    secs = floor(Int, s % 60)
    mins = floor(Int, s / 60.0)

    h = floor(Int,mins / 60)

    mins = mins % 60

    string("Time: ", h, ":", mins, ":", secs)

end

function simulationRound(body, systemMatrix, τ, h2, α)

    result =  (τ / h2) * (α .* systemMatrix) * body  .+ body

    result
end


function simulation()
    α = 3.8e-6 # m^2 / s
    temp_danger = 45.0 # °C
    length = 1.0 # m
    measure_danger = (80) # idx

    grid_size = (100) 

    flame_temp = 100.0 # °C
    flame_location = 1 # idx


    #τ = 0.01 # ?
    h = length / grid_size[1]
    h2 = h * h
    τ = h2 / (2 * α)

    ##

    stick = zeros(grid_size)

    stick[flame_location] =  flame_temp;
    t = 0
    m = bandMatrix(grid_size[1])

    
    #while stick[measure_danger[1]] < temp_danger
        stick .= simulationRound(stick, m, τ, h2,  α)
        (t % 20000 == 0) && print(t,τ, stick[measure_danger[1]])

        t += 1
    #end

    println(secondsToTime(t * τ))
end

function print(t, τ, temp) 
    println(secondsToTime(t * τ), " Temp: " , temp)
    false
end

function bandMatrix(size)
    m = zeros((size,size));
    #m[1,1:2] .= [-2 1]
    for i in 2:size-1
        m[i, i-1:i+1] = [1 -2 1]
    end
    m[size, end-1:end] = [1 -1]
    m
end


#@show Threads.nthreads()
simulation()
