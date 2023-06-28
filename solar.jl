using Printf

B_GregBias(Y::Int; isGreg=true) = !isGreg ? 0 : (2 - (Y÷100) + (Y÷400))
monthconvert(Y::Int, M::Int) = M > 2 ? (Y,M) : (Y-1, M+12)
C_J = 36525
Δ_JD0 = 2451545
C_T = 1.002738
C_T0 = 2400.05134
C_Θ = 6.697376
JulianDate(Y::Int,M::Int,D::Int; isGreg=true) = floor((C_J/100)*(Y + 4716)) + floor(30.6001(M+1)) + D + B_GregBias(Y, isGreg=isGreg) - 1524.5

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
azimuth(δ, ϕ, τ) = atan(sin(τ)/ (cos(τ) * sin(ϕ) - tan(δ) * cos(ϕ)))
elevation(δ, ϕ ,τ) = asin(cos(δ) * cos(τ) * cos(ϕ) + sin(δ) * sin(ϕ))



function calc(Y::Int, M::Int, D::Int, T::Float64) 
    long = 11.6 # Geographische laenge
    lat = deg2rad(48.1) # Geografische Breite
    δ = 0.0  # TODO Deklination der Sonne; 
    α = 0.0 # TODO Right Ascention ?= Rektatention

    (Y,M) = monthconvert(Y,M)
    jd0 = JulianDate(Y,M,D)
    @show jd0 
    temp = (jd0 - Δ_JD0)
    n = temp + T 
    t0 = temp / C_J
    @show t0 typeof(t0)
    #t0 = 0.06594113621
    Θ_hG = C_Θ + C_T0 * t0 + C_T * T * 24.0
    Θ_hG  %= 24
    @show Θ_hG
    Θ_G = Θ_hG * 15
    Θ = Θ_G + long
    @show Θ
    
    
    (α, δ) = SunProjection(n)
    #(α, δ)
    τ = deg2rad(Θ) - α 
    #(τ, δ)

    (α, δ, azimuth(δ,lat,τ), elevation(δ,lat,τ))
end


@printf "Rektaszension: %.4f, Deklination: %.4f\nazimuth: %.2f, elevation: %.2f\n" (calc(2006,8,6,0.25) .|> rad2deg)...
# @printf "azimuth: %.2f, elevation: %.2f\n" calc(2023,6,28,0.5)...
