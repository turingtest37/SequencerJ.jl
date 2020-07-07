using Colonies
using SequencerJulia: kl, EMD
using Distances: pairwise

function dbb(m,n,f)
    A = m .* n
    @debug "A" A
    B = pairwise(EMD(), A; dims=1)
    # B = kl(A)
    @debug "B" B
    kooky(i) = oftype(i,rand(-10:-1))
    C = map(i->isinf(i) || isnan(i) ? 0 : i, B)
    @debug "C" C
    sum(C)
    # "sum(map(i->isinf(i) ? sign(i) * oftype(i,1) : i, kl(m .* n)))"
end

generatemany(10,10,50,50,true;
       limit = 20,
       reducef = ColonyFunction("sum(map(i->isinf(i) ? 0 : i, pairwise(EMD(), m .* n; dims=1)))", dbb),
       maskdim = 3,
       seed = seedwith(square),
       destdir = "/Users/doug/dev/automaton/img")
