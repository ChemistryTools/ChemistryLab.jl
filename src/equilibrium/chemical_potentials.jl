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
        _n = max.(n, p.ϵ)
        ΔₐG⁰overT = similar(_n)
        ΔₐG⁰overT .= p.ΔₐG⁰overT
        ΔₐG⁰overT[idxsolute] .+= log(Csolvent)
        ntot_aqueous = _n[idxsolvent] + sum(_n[idxsolute])
        if !iszero(ntot_aqueous)
            ΔₐG⁰overT[idxsolute] .+= log.(_n[idxsolute] ./ _n[idxsolvent])
            ΔₐG⁰overT[idxsolvent] += log(_n[idxsolvent] / ntot_aqueous)
        end
        ntot_gas = sum(_n[idxgas])
        if !iszero(ntot_gas)
            ΔₐG⁰overT[idxgas] .+= log.(_n[idxgas] ./ ntot_gas)
        end
        return ΔₐG⁰overT
    end
    return μ
end
