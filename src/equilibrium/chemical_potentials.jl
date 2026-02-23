@kwdef struct State{S<:AbstractSpecies,TS<:Number,Ttype<:Number,Ptype<:Number,ttype<:Number}
    species::OrderedDict{S,TS} = OrderedDict{AbstractSpecies, Number}()
    T::Ttype = 298.15u"K"
    P::Ptype = 1.0u"bar"
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

function potentials_dilute_ideal(species::AbstractVector)
    idxsolvent = findfirst(x->x.class == SC_AQSOLVENT, species)
    @assert !isnothing(idxsolvent) "No solvent has been found in the mix"
    idxsolute = findall(x -> x.class == SC_AQSOLUTE, species)
    idxcrystal = findall(x -> x.aggregate_state == AS_CRYSTAL, species)
    idxgas = findall(x -> x.aggregate_state == AS_GAS, species)
    R = ustrip(Constants.R)
    Csolvent = ustrip(1/species[idxsolvent].M)# 55.5 only for water solvent
    function μ(n, p)
        p = NamedTuple(p)
        T = ustrip(get(p, :T, 298.15))
        RT = R * T
        _n = max.(n, p.ϵ)
        ΔₐGpoverRT = similar(_n)
        ΔₐGpoverRT .= p.ΔₐG⁰/RT
        ΔₐGpoverRT[idxsolute] .+= log(Csolvent)
        ntot_aqueous = _n[idxsolvent] + sum(_n[idxsolute])
        if !iszero(ntot_aqueous)
            ΔₐGpoverRT[idxsolute] .+= log.(_n[idxsolute] ./ _n[idxsolvent])
            ΔₐGpoverRT[idxsolvent] += log(_n[idxsolvent] / ntot_aqueous)
        end
        ntot_gas = sum(_n[idxgas])
        if !iszero(ntot_gas)
            ΔₐGpoverRT[idxgas] .+= log.(_n[idxgas] ./ ntot_gas)
        end
        return ΔₐGpoverRT
    end
    return μ
end
