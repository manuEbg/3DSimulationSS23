B_GregBias(Y::Int; isGreg=true) = !isGreg ? 0 : (2 - (Y÷100) + (Y÷400))
M_c(M::Int, Y::Int) = M > 2 ? (Y,M) : (Y-1, M+12)
C_J = 36525
Δ_JD0 = 2451545
C_T = 1.002738
C_T0 = 2400.05134
C_Θ = 6.697376
JD_0(Y::Int,M::Int,D::Int; isGreg=true) = floor((C_J/100)(Y + 4716)) + floor(30.6001(M+1)) + D + B_GregBias(Y, isGreg=isGreg) - 1524.5
T_0(JD_0) = (JD_0 - Δ_JD0)/C_J
ΘhG(Y::Int, M::Int, D::Int, T; isGreg=true) = CΘ + C_T0 * T_0(JD_0(M_c(M, Y)..., D, isGreg=isGreg)) + C_T * T
Θ_hG(2023, 6, 14, 18.34)