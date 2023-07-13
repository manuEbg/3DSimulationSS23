using Plots
using Printf
using LinearAlgebra
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

# TODO support something other than UT (Universal Time)
struct Date
  T::Float64 # Fraction of Day 0.0..1.0
  T_H::Float64 # Fraction of Day 0.0..24.0
  MJD::Float64
  MJD0::Float64
  """
    t : Fraction of Day. 
    Note: 0.0 => 00:00, 0.5 => 12:00, 1.0 => 00:00 on the next day
  """
  function Date(y::Int, m::Int, d::Int, t::Float64) 
    jd0 = juliandate(y,m,d)
    mjd0 = jd0 - Δ_JD0
    mjd = mjd0 + t # n

    new(t,t*24.0,mjd,mjd0)
  end

  function Date(y::Int, m::Int, d::Int)
    new(y,m,d,0.0)
  end

  """
    Date(date::Date, Δt::Number)

    constructs a new date from the given date and adds Δt (a value of 1.0 results in the next day)
  """
  Date(date::Date, Δt::Number) = new(date.T + Δt, date.T_H + Δt * 24, date.MJD + Δt, date.MJD0)

  """
    JulianDate(Y::Integer, M::Integer, D::Integer; isGreg=true)
    Returns the JulianDate on the given day at UT 00:00
  """
  function juliandate(Y::Integer, M::Integer, D::Integer; isGreg=true)::Float64
      (y, m) = monthconvert(Y, M)
      bias = isGreg ? B_GregBias(y) : 0
      century = floor((C_J / 100) * (y + 4716))
      mon = floor(30.6001(m + 1))
      century + mon + D + bias - 1524.5
  end
end

Base.:+(a::Date, b::Number) = Date(a, b)

struct GeoLocation 
  longitude::Angle
  latitude::Angle
end

struct PanelSize
  width::Float64
  height::Float64
end

struct SphericalCoordinates
  azimuth::Angle
  elevation::Angle
  SphericalCoordinates() = new(0 |> deg, 0 |> deg)
  SphericalCoordinates(azimuth, elevation) = new(azimuth,elevation)
end

Base.:-(a::SphericalCoordinates, b::SphericalCoordinates) = SphericalCoordinates(
  a.azimuth - b.azimuth, 
  a.elevation - b.elevation 
)

function to_cartesian(pos::SphericalCoordinates)
  ϕ = (pos.azimuth + (180 |> deg)).rad
  Θ = (π/2) - pos.elevation.rad
  [sin(Θ) * cos(ϕ), sin(Θ) * sin(ϕ), cos(Θ)]
end

"""
  SunProjection 
  α - Right Ascension ( The angle of the sun projected onto the equatorial plane. 
      Measured in counter clockwise direction from the vernal equinox )
  δ - Declination ( The Latitude of the projection of the sun )
"""
struct SunProjection
  α::Angle 
  δ::Angle 
  
  function SunProjection(N)
    Ε = 23.439 - 0.4e-6 * N |> deg2rad # angle between earth rotation axis and eclipsis
    g = 357.528 + 0.9856003 * N |> deg2rad
    L = 280.460 + 0.985647 * N # degree
    Λ = L + 1.915 * sin(g) + 0.01997 * sin(2 * g) |> deg2rad
    α = atan(cos(Ε) * tan(Λ)) # Right Ascension
  
    α = cos(Λ) <= 0 ? α + 4 * atan(1) : α
  
    δ = asin(sin(Λ) * sin(Ε)) # Declination
  
    new(α |> rad, δ |> rad)
  end
end

"""
All values as radiants
-  δ Deklination der Sonne
-  ϕ Geografische Breite Zielort
-  τ Winkeldistanz Zielort Sonne
"""
azimuth(δ, ϕ, τ) = atan(sin(τ) , cos(τ) * sin(ϕ) - tan(δ) * cos(ϕ))
elevation(δ, ϕ, τ) = asin(cos(δ) * cos(τ) * cos(ϕ) + sin(δ) * sin(ϕ))

function sun_spherical(date::Date, location::GeoLocation; refraction_correction::Function = identity)::SphericalCoordinates
  sun = SunProjection(date.MJD)
  τ = Θ(date,location) - sun.α
  SphericalCoordinates(
    azimuth(sun.δ.rad, location.latitude.rad, τ.rad) |> rad,
    elevation(sun.δ.rad, location.latitude.rad, τ.rad) |> rad |> refraction_correction
  )
end

struct SolarPanel
  location::GeoLocation
  position::SphericalCoordinates
  size::PanelSize
end  

function Θ(date::Date, location::GeoLocation)::Angle
  t0 = date.MJD0 / C_J
  Θ_hG = C_Θ + C_T0 * t0 + C_T * date.T_H
  Θ_hG %= 24
  Θ_G = deg(Θ_hG * 15)
  Θ_G + location.longitude
end

function power_factor(sun:: SphericalCoordinates, panel:: SphericalCoordinates)::Float64 
  f = dot(to_cartesian(sun), to_cartesian(panel))
  f = (sun.elevation.deg > 0 && f >= 0) * f
  f
end

"""
    refractionCorrection(elevation::Angle)::Angle
    - elevation of the sun
    - returns corrected elevation in degree
    due to the refraction of the atmosphere the suns elevation needs to be adjusted.
"""
function refractionCorrection(elevation::Angle)::Angle
  h = elevation.deg
  R = 1.02 / tand(h + 10.3 / (h + 5.11))
  deg(h + R / 60)
end

struct Simulation
  solarpanel::SolarPanel
  startdate::Date
  enddate::Date
  timestep::Float64 # timestep in days
end

struct SimulationState
  t::Date # Current simulation time
  #sun_eq::SunProjection # 
  sun_sph::SphericalCoordinates
  sun_power::Float64
  panel_power::Float64
  power_consumption::Float64
  storage_energy::Float64
  SimulationState(date::Date) = new(date, SphericalCoordinates(), 0,0,0,0)
  SimulationState(t,sun_sph,sun_power,panel_power,power_consumption,storage_energy) = new(
    t,sun_sph,sun_power,panel_power,power_consumption,storage_energy
  )
end


"""
    (sim::Simulation)(state::SimulationState)::SimulationState

One step forward in the simulation
"""
function (sim::Simulation)(state::SimulationState)::SimulationState
    date = state.t + sim.timestep
    sun = sun_spherical(
      date, 
      sim.solarpanel.location, 
      refraction_correction = refractionCorrection
    )

    # sunrise and sunset
    -0.1 < sun.elevation.deg < 0.1 && @show date.T_H 

    sun_p = 1000 # Watts per Square Meter
    panel_p = power_factor(sun, panel.position) * 205 # Watts

    e = state.storage_energy + panel_p * sim.timestep # Watts/Day ?
    consumption = state.power_consumption * sim.timestep
    e = e - consumption

    SimulationState(date,sun,sun_p,panel_p, consumption, e)
end

"""
    (sim::Simulation)(data_collection_freq = 10)::Array{T,2}

  run the simulation.

  collect plotting date every `data_collection_freq` steps
"""
function (sim::Simulation)(data_collection_freq = 10)::Array{Float64,2}
  data = Array{Float64}(undef, 4, 0)
  state = SimulationState(sim.startdate)
  i = 1
  while(state.t.MJD < sim.enddate.MJD)
    state = sim(state)
    if i % data_collection_freq == 0
      data = hcat(data, [
      state.sun_sph.azimuth.deg, 
      state.sun_sph.elevation.deg, 
      state.panel_power, 
      state.storage_energy]
      )
    end
    i = i + 1
  end
  data
end

function visualize(data::Array{Float64, 2})
  plot(
    1:size(data,2), 
    data', 
    lab = ["SunAz" "SunH" "Power" "Energy"], 
    legend = :outerright
    )
end

panel = SolarPanel(
  GeoLocation(11.6 |> deg, 48.1 |> deg), # Munich
  SphericalCoordinates(57.34 |> deg , 90 |> deg), 
  PanelSize(1,1)
)

simulation = Simulation(
  panel,
  Date(2023,6,20),
  Date(2023,6,23), 
  0.01
)

data = simulation(10)
@printf "simulation finished\n"

pl = visualize(data)
Plots.pdf(pl, "testr.pdf")
@printf "Plot printed\n"



