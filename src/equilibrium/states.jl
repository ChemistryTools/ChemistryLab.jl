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
