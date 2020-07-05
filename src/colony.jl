using Colonies
using SequencerJulia

dbb(m,n,f) = sum(m .* map(i->isinf(i) ? Int(1) : Int(i), kl(n)))

generatemany(10,10,50,50,true;
       limit=10,
       reducef = ColonyFunction(dbb),
       maskdim=5,
       destdir="/home/doug/dev/automaton/img")
