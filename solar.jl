using Printf
using Plots

B_GregBias(Y::Int) = (2 - (Y÷100) + (Y÷400))
monthconvert(Y::Int, M::Int) = M > 2 ? (Y,M) : (Y-1, M+12)
C_J = 36525
Δ_JD0 = 2451545
C_T = 1.002738
C_T0 = 2400.05134
C_Θ = 6.697376
function JulianDate(Y::Int, M::Int, D::Int; isGreg=true) 
    (y,m) = monthconvert(Y,M)
    bias = isGreg ? B_GregBias(y) : 0
    century = floor((C_J/100)*(y + 4716))
    mon = floor(30.6001(m+1))
    century + mon + D + bias  - 1524.5
end
function SunProjection(N)
    Ε = 23.439 - 0.4e-6 * N |> deg2rad
    g = 357.528 + 0.9856003 * N |> deg2rad
    L = 280.460 + 0.985647 * N # degree
    Λ = L + 1.915 * sin(g) + 0.01997 * sin(2*g) |> deg2rad
    α = atan(cos(Ε) * tan(Λ)) # Rektaszension

    if (cos(Λ) <= 0)
        α = α + 4 * atan(1)
    end

    δ = asin(sin(Λ) * sin(Ε)) # Deklination

    (α,δ)
end




# δ Deklination der Sonne
# ϕ Geografische Breite Zielort
# τ Winkeldistanz Zielort Sonne
sunPosition(δ, ϕ, τ) = (azimuth(δ,ϕ,τ), elevation(δ,ϕ,τ))
azimuth(δ, ϕ, τ) = atan(sin(τ)/ (cos(τ) * sin(ϕ) - tan(δ) * cos(ϕ)))
elevation(δ, ϕ ,τ) = asin(cos(δ) * cos(τ) * cos(ϕ) + sin(δ) * sin(ϕ))

# longitude of the destination in radiants
function calc_n_and_Θ(jd0, longitude, T)
    temp = (jd0 - Δ_JD0)
    n = temp + T 
    t0 = temp / C_J
    Θ_hG = C_Θ + C_T0 * t0 + C_T * T * 24.0
    Θ_hG  %= 24
    Θ_G = Θ_hG * 15 |> deg2rad
    
    (n, Θ_G + longitude)
end

function calc(Y::Int, M::Int, D::Int, T::Float64) 
    long = 11.6 |> deg2rad # Geographische laenge
    lat = 48.1 |> deg2rad # Geografische Breite

    jd0 = JulianDate(Y,M,D)
    
    (n, Θ) = calc_n_and_Θ(jd0, long, T)

    (α, δ) = SunProjection(n)

    τ = Θ - α 
    #@show Θ |> rad2deg

    (a,h) = sunPosition(δ,lat,τ)
    #@show h |> rad2deg
    h = refractionCorrection(h)
    (a, h)
end

# elevation as radiants
function refractionCorrection(elevation::Float64)
    h = elevation |> rad2deg
    R = 1.02 / tand(h + 10.3/(h + 5.11))
    (h + R/60) |> deg2rad
end
range = 0:0.05:800.0
result = [calc(2006,8,6,t/24.0) .|> rad2deg for t in range]

function printRes(result)
    for r in result
        @printf "azimuth: %.4f, elevation: %.4f\n" r...
    end
end


data = (vec([r[2] for r in result]), vec([r[1] for r in result])) |> a -> hcat(a...)

    
plot(range, data[:,1])
# Rektaszension: %.4f, Deklination: %.4f\n