@kwdef struct State{S<:AbstractSpecies,TS<:Number,Ttype<:Number,Ptype<:Number,ttype<:Number}
    species::OrderedDict{S,TS} = OrderedDict{AbstractSpecies, Number}()
    T::Ttype = 298.15u"K"
    P::Ptype = 1.0u"bar"
    t::ttype = 0.0u"s"
end

function Base.show(io::IO, st::State)
    println(io, "State at T: $(st.T) ◆ P: $(safe_uconvert(us"bar", st.P)) ◆ t: $(st.t)")
    print(st.species)
end

function Base.show(io::IO, ::MIME"text/plain", st::State)
    println(io, "State at T: $(st.T) ◆ P: $(safe_uconvert(us"bar", st.P)) ◆ t: $(st.t)")
    pad = max(length.(symbol.(keys(st.species)))...)
    for (k,v) in st.species
        println(io, lpad(symbol(k), pad), ": ", v)
    end
end

function potentials_dilute(species::AbstractVector{<:AbstractSpecies})
    idxsolvent = findfirst(x->x.class == SC_AQSOLVENT, species)
    @assert !isnothing(idxsolvent) "No solvent has been found in the mix"
    R = ustrip(Constants.R)
    function μ(n, p)
        pp = NamedTuple(p)
        T = ustrip(pp.T)
        RT = R * T
        _n = max.(n, pp.ϵ)
        # 55.5 only for water solvent
        ΔₐGpoverRT = pp.ΔₐG⁰/RT + log(55.5)*[i == idxsolvent ? 0 : 1 for i in eachindex(n)]
        ntot = sum(_n)
        return ΔₐGpoverRT + log.(_n ./ ntot)
    end
    return μ
end
