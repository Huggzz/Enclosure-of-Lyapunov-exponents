using Polynomials, PolynomialRoots, IntervalArithmetic, Combinatorics, Serialization

prec = 2048
setprecision(prec)

M = big(550)

z = Polynomial(interval.(BigFloat, [0, 1]))
Hprev = Polynomial([interval(BigFloat, 1)])
H = Polynomial(interval.(BigFloat, [0, 2]))

sups = zeros(Interval{Float64}, M+2)

sups[1] = interval(BigFloat, 1.0)

ϵ  = 2^(-prec//2)

for m = 1:M+1
    println(m)
    DH = interval(BigFloat, 2*m)*Hprev - z*H
    X = sort(real.(PolynomialRoots.roots(mid.(DH.coeffs))))
    Xl = interval.(X .- ϵ)
    Xu = interval.(X .+ ϵ)
    if all(sign.(DH.(Xl)).*sign.(DH.(Xu)) .== -1) && all(sup.(Xl[1:end-1]).< inf.(Xu[2:end])) && length(X) == m + 1
        Xrig = interval.(inf.(Xl), sup.(Xu))
        sups[m+1] = interval(Float64, interval(BigFloat,maximum(sup.(abs.(H.(Xrig)).*exp.(-Xrig.^2/interval(2)))))/sqrt(interval(BigFloat, 2^m*factorial(m))))
    else
        println("problem")
    end
    H, Hprev = interval(2)*z*H - interval(BigFloat, 2*m)*Hprev, H
end

serialize("pendulum_sups", sups)