export forward1, central2

function forward1(x::Vector, h::Real)
    return [(x[2:end] - x[1:end-1]) / h..., (x[end] - x[end-1]) / h]
end

function forward1(x::Vector, t::Vector)
    return [((x[2:end] - x[1:end-1]) ./ (t[2:end] - t[1:end-1]))..., (x[end] - x[end-1]) / (t[end] - t[end-1])]
end

function central2(x::Vector, h::Real)
    return [(x[2] - x[1]) / h, (x[3:end] - x[1:end-2]) / 2h..., (x[end] - x[end-1]) / h]
end

function central2(x::Vector, t::Vector)
    return [(x[2] - x[1]) / (t[2] - t[1]), ((x[3:end] - x[1:end-2]) ./ (t[3:end] - t[1:end-2]))..., (x[end] - x[end-1]) / (t[end] - t[end-1])]
end