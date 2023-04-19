## Helpers

literPerMin_to_cm3PerS(x) = (x * 1000.0) / 60.0
literPerMin_to_m3PerS(x) = (x / 60000.0)

## Constants

g = 9.81 # m/s^2

## the Container

mutable struct Container
    length::Float64
    width::Float64
    height::Float64
    area::Float64
    outputDistance::Float64
    nOutputs::Int
    outputHeight::Vector{Float64}
    outputArea::Vector{Float64}
    waterLevel::Float64


    function Container(l,w,h,outputDistance::Float64) 

        ## Calculate the Area of the outputs

        # 5.32 mm^2 ist die flaeche der loecher
        targetSpeed = literPerMin_to_m3PerS(1.0)
        waterColumn = 0.5 # m
        outputArea = targetSpeed / sqrt(waterColumn * g * 2) # m

        n = floor(Int, h / outputDistance) + 1
        oh = vec([i * outputDistance for i in 0:n-1])
        oa = vec([outputArea for i in 1:n])
        new(
            l,w,h,
            l*w,
            outputDistance,
            n,
            oh,
            oa,
            0.0 
        )
    end
end

function nActiveOutputs(c::Container)
    floor(Int, c.waterLevel / c.outputDistance) + 1
end

function increaseWaterLevel!(c::Container, ϕ::Float64, t::Float64)
    c.waterLevel = c.waterLevel + t * ϕ / c.area
end

function decreaseWaterLevel!(c::Container, ϕ::Float64, t::Float64) 
    
    h = c.waterLevel - ϕ * t / c.area

    c.waterLevel = h > 0 ? h : 0
end

function getOutputSpeed(c::Container)
    ϕ = 0.0;
    for i in 1:c.nOutputs
        waterCol = c.waterLevel - c.outputHeight[i]

        speed = waterCol > 0.0 ? c.outputArea[i] * sqrt(waterCol * g * 2) : 0.0

        ϕ = ϕ + speed
    end
    ϕ
end

## Simulation
function simulate(c::Container, ϕ_in::Float64, timestep::Float64)
    time = 0
    while true
        
        increaseWaterLevel!(c, ϕ_in, timestep)

        if(c.waterLevel >= 0.9) 
            break; 
        end
        ϕ_out = getOutputSpeed(c)
        
        decreaseWaterLevel!(c, ϕ_out, timestep)

        time = time + 1;
        if(time % 30000 == 0) 
            level = c.waterLevel
            t = time * timestep
            @show t level
        end
        
    end


    c.waterLevel = 0

    @show time * timestep
end



## Configuration
ϕ_in = literPerMin_to_m3PerS(10)
#ϕ_out = literPerMin_to_m3PerS(1)
c = Container(0.1,0.1,1.0,0.1)
timestep = 0.001;
simulate(c, ϕ_in, timestep)

