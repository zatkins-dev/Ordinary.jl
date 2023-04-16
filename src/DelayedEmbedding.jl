export dce

function myrem(x, y)
    x - y * round(x / y)
end

function dce(x, h, m, τ)
    if (myrem(τ, h) > eps(Float64))
        error("τ should be a multiple of h")
    end
    k = round(Int, τ / h)
    reduce(hcat, [[x[j+i*k] for i in 0:m-1] for j in 1:length(x)-k*(m-1)-1])
end