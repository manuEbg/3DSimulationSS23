using LinearAlgebra
using PlotlyJS;

function secondsToTime(s)

    secs = floor(Int, s % 60)
    mins = floor(Int, s / 60.0)

    h = floor(Int,mins / 60)

    mins = mins % 60

    string("Time: ", h, ":", mins, ":", secs)

end

function simulationRound(body, systemMatrix, border_cases, τ, h2, α)

    result =  (τ / h2) * ((α .* systemMatrix) * body)  .+ body .+ border_cases

    result 
end


function simulation()
    
    α = 3.8e-6 # m^2 / s
    temp_danger = 45.0 # °C
    length = 1.0 # m
    
    grid_size = (5, 50) 
    measure_danger = floor(Int,(0.8 * grid_size[2]) + grid_size[2] * floor(Int,grid_size[1]/2) ) # idx
    @show measure_danger
    flame_temp = 100.0 # °C
    flame_location = 1 # idx

    h = length / grid_size[2]
    h2 = h * h
    τ = h2 / (2 * α)
    @show τ
    ##

    stick = zeros(grid_size)

    stick[flame_location, :] .=  flame_temp;
    t = 0
    m = generateSystemMatrix(stick)
    
    #

    #return m
    v_stick = vec(stick);
    
    border_cases = zeros(size(v_stick))
    #while v_stick[measure_danger] < temp_danger
    for avc in 1:1
        v_stick .= simulationRound(v_stick, m, border_cases, τ, h2,  α)
        (t % 2 == 0) && print(t,τ, v_stick[measure_danger])
        t += 1
    end
    
    println(secondsToTime(t * τ))
    stick = reshape(v_stick, grid_size)
    stick
end

function print(t, τ, temp) 
    println(secondsToTime(t * τ), " Temp: " , temp)
    false
end
 
function generateSystemMatrix(grid)

    # grid = reshape([i for i in 1:12],(3, 4))

    v = vec(grid)

    sysMat = zeros((length(v), length(v)))


    row_count = size(grid,1);
    col_count = size(grid,2);

    for (i,node) in pairs(v)
        row = ((i-1) % row_count) + 1
        col = ceil(Int, i / float(row_count))


        if row == 1 
            sysMat[i, i:i+1] = [-4 1]
        elseif row == row_count
            sysMat[i, i-1:i] = [1 -4]   
        else
            sysMat[i, i-1:i+1] = [1 -4 1]
        end

        if col == col_count
            sysMat[i, i-row_count] = 1 
        elseif col == 1
            sysMat[i, i+row_count] = 1   
        else
            sysMat[i, i+row_count] = 1   
            sysMat[i, i-row_count] = 1
        end  
        
    end
    dirichlet!(sysMat, [1])
    sysMat
end


function dirichlet!(mat, indices)
    mat[indices, :] .= 0
end

##



#plot(heatmap(z=sysMat))


##

#@show Threads.nthreads()
s = simulation()
@show s
plot(heatmap(z=s))