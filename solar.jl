using Plots
using Printf
# https://static.trinasolar.com/sites/default/files/DE_Datasheet_VertexS_DE09.08_2021_A.pdf
B_GregBias(Y::Int) = (2 - (Y ÷ 100) + (Y ÷ 400))
monthconvert(Y::Int, M::Int) = M > 2 ? (Y, M) : (Y - 1, M + 12)
C_J = 36525
Δ_JD0 = 2451545 # JD on 01.01.2000 00:00 (J2000)
C_T = 1.002738 # Factor for fraction of day
C_T0 = 2400.05134 #
C_Θ = 6.697376 # Startime in greenwich on J2000

struct Angle 
  rad
  deg
end
deg(value) = Angle(value |> deg2rad, value)
rad(value) = Angle(value, value |> rad2deg)
Base.:+(a::Angle,b::Angle) = rad(a.rad + b.rad)
Base.:-(a::Angle,b::Angle) = rad(a.rad - b.rad)
struct Date
  Y::Int
  M::Int
  D::Int
  T::Float64
  T_H::Float64
  MJD::Float64
  MJD0::Float64
  """
    t : Fraction of Day. 
    Note: 0.0 => 00:00, 0.5 => 12:00, 1.0 => 00:00 on the next day
  """
  function Date(y::Int, m::Int, d::Int, t::Float64) 
    jd0 = JulianDate(y,m,d)
    mjd0 = jd0 - Δ_JD0
    mjd = mjd0 + t # n

    new(y,m,d,t,t*24.0,mjd,mjd0)
  end

  function JulianDate(Y::Int, M::Int, D::Int; isGreg=true)
    (y, m) = monthconvert(Y, M)
    bias = isGreg ? B_GregBias(y) : 0
    century = floor((C_J / 100) * (y + 4716))
    mon = floor(30.6001(m + 1))
    century + mon + D + bias - 1524.5
  end

end

struct GeoLocation 
  longitude::Angle
  latitude::Angle
end

struct PanelSize
  width::Float64
  height::Float64
end

struct AnglePosition
  azimuth::Angle
  elevation::Angle
end

struct SunProjection
  α::Angle # Rektaszension
  δ::Angle # Deklination
  function SunProjection(N)
    Ε = 23.439 - 0.4e-6 * N |> deg2rad # angle between earth rotation axis and eclipsis
    g = 357.528 + 0.9856003 * N |> deg2rad
    L = 280.460 + 0.985647 * N # degree
    Λ = L + 1.915 * sin(g) + 0.01997 * sin(2 * g) |> deg2rad
    α = atan(cos(Ε) * tan(Λ)) # Rektaszension
  
    α = cos(Λ) <= 0 ? α + 4 * atan(1) : α
  
    δ = asin(sin(Λ) * sin(Ε)) # Deklination
  
    new(α |> rad, δ |> rad)
  end
end

"""
All valuse as radiants
-  δ Deklination der Sonne
-  ϕ Geografische Breite Zielort
-  τ Winkeldistanz Zielort Sonne
"""
azimuth(δ, ϕ, τ) = atan(sin(τ) , cos(τ) * sin(ϕ) - tan(δ) * cos(ϕ))
elevation(δ, ϕ, τ) = asin(cos(δ) * cos(τ) * cos(ϕ) + sin(δ) * sin(ϕ))

function sun_angles(date::Date, location::GeoLocation; refraction_correction::Function = identity)::AnglePosition
  sun = SunProjection(date.MJD)
  τ = Θ(date,location) - sun.α
  AnglePosition(
    azimuth(sun.δ.rad, location.latitude.rad, τ.rad) |> rad,
    elevation(sun.δ.rad, location.latitude.rad, τ.rad) |> refraction_correction |> rad
  )
end
struct SolarPanel
  location::GeoLocation
  position::AnglePosition
  size::PanelSize
end  

function Θ(date::Date, location::GeoLocation)::Angle
  t0 = date.MJD0 / C_J
  Θ_hG = C_Θ + C_T0 * t0 + C_T * date.T_H
  Θ_hG %= 24
  Θ_G = deg(Θ_hG * 15)
  Θ_G + location.longitude
end

function calc(date::Date)
  l = GeoLocation(11.6 |> deg, 48.1 |> deg)
  
  sun = sun_angles(date, l, refraction_correction = identity)
  # h = refractionCorrection(h)
  sun
end

# elevation as radiants
"""
    refractionCorrection(elevation::Float64)::Float64
    - elevation in degree
    - returns corrected elevation in degree
    due to the refraction of the atmosphere the suns elevation needs to be adjusted.
"""
function refractionCorrection(elevation::Float64)::Float64
  h = elevation
  R = 1.02 / tand(h + 10.3 / (h + 5.11))
  (h + R / 60)
end
range = 0:0.05:72.0

data = [calc(Date(2023,9,21,t/24.0)) for t in range]
plotting_data = hcat(
  map(sun -> sun.azimuth.deg , data),
  map(sun -> sun.elevation.deg, data) 
)

function printRes(result)
  for r in result
    @printf "azimuth: %.4f, elevation: %.4f\n" r.azimuth.deg r.elevation.deg
  end
end

# printRes(data)

pl = plot(range, plotting_data)
Plots.pdf(pl, "test.pdf")
