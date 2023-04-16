
unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))