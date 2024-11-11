using IntervalArithmetic, Combinatorics, Polynomials, PolynomialRoots, Serialization, Random, BlockArrays, Base.Threads

prec = 2048
setprecision(prec)
Ns = collect(500:-2:0)
NS = vcat(Ns[1]+1, Ns);
rows = reverse(sort(vcat(NS[1]+3, NS.+1, NS.+1)))

sups = BlockVector(zeros(Interval{BigFloat}, (sum(rows))), rows)
ϵ = interval(2^(-prec//2))
println("computing sups")


Threads.@threads for k in shuffle(big.(collect(1:length(NS))))
    supsk = zeros(Interval{BigFloat}, length(sups[Block(2*k)]))
    println((Threads.threadid(),k))
    for m=0:big(NS[k])
        # println((k,m))
        Lₘ = Polynomial([interval(-1)^interval(j)*interval.(BigFloat, binomial(m+2*k,m-j)//factorial(j)) for j=0:m])
        if m>0
            lₘ = Polynomial([interval(-1)^interval(j)*interval.(BigFloat, binomial(m-1+2*k+1,m-1-j)//factorial(j)) for j=0:m-1])
        else
            lₘ = Polynomial([interval(BigFloat,0)])
        end
        P = Polynomial([interval(0), interval(-1)])*lₘ+Polynomial([interval(k), interval(0), interval(-1)])*Lₘ
        x = interval.(sort(real.(PolynomialRoots.roots(mid.(P.coeffs)))))
        xl = x .- ϵ
        xu = x .+ ϵ

        if all(sign.(P.(xl)).*sign.(P.(xu)) .==-1) && length(x) == m+2 && all(sup.(xu[1:end-1]) .< inf.(xl[2:end]))
            Z = sqrt(interval.(BigFloat, factorial(m+2*k)//factorial(m)))
            xrig = interval.(inf.(xl[sup.(xu) .>= 0]), sup.(xu[sup.(xu) .>= 0]));
            supsk[m+1] = interval(maximum(sup.(abs.(Lₘ.(xrig).*xrig.^k.*exp.(-xrig.^2/interval(2))))))/Z
        else
            println((k, m))
        end
    end
    sups[Block(2*k)] = sups[Block(2*k+1)] = supsk
end


for m=0:big(NS[1]+2)
    # println((k,m))
    Lₘ = Polynomial([interval(-1)^interval(j)*interval.(BigFloat, binomial(m,m-j)//factorial(j)) for j=0:m])
    if m>0
        lₘ = Polynomial([interval(-1)^interval(j)*interval.(BigFloat, binomial(m-1+1,m-1-j)//factorial(j)) for j=0:m-1])
    else
        lₘ = Polynomial([interval(BigFloat,0)])
    end
    P = lₘ+Polynomial([interval(0), interval(1)])*Lₘ
    x = interval.(sort(real.(PolynomialRoots.roots(mid.(P.coeffs)))))
    xl = x .- ϵ
    xu = x .+ ϵ

    if all(sign.(P.(xl)).*sign.(P.(xu)) .==-1) && length(x) == m+1 && all(sup.(xu[1:end-1]) .< inf.(xl[2:end]))
        xrig = interval.(inf.(xl[sup.(xu) .>= 0]), sup.(xu[sup.(xu) .>= 0]));
        sups[Int64(m+1)] = interval(maximum(sup.(abs.(Lₘ.(xrig).*exp.(-xrig.^2/interval(2))))))
    else
        println(m)
    end
end

Threads.@threads for k in shuffle(big.(collect(1:length(NS))))
    sups[Block(2*k)] = sups[Block(2*k+1)]
end

serialize("suppsi", interval.(Float64, sups))