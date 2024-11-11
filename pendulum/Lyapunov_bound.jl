using IntervalArithmetic, Serialization

prec = 256
setprecision(prec)

γ = interval(2//3)
κ = interval(Float64, 1//4)
σ = interval(4)
c = interval(Float64, 185//1000)

ylim = (interval(γ)+sqrt(interval(γ)^2+(interval(2)*c+κ)*σ^2))/κ

Y = mince(interval(0, sup(ylim)), 1000000)

d = interval(maximum(sup.((κ*(-κ*Y.^2 .+σ^2+interval(2)*γ*Y)/(interval(2)*σ^2).+c).*exp.(κ*Y.^2/(2*σ^2)))))

serialize("pendulum_Lyapunov_bound", d/c)
