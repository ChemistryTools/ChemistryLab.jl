abstract type Callable end  # Base type for all callable thermodynamic functions

struct BasisFunction{D,I} <: Callable end

const ExtendedDegree = Union{Int,Symbol}
degree(::BasisFunction{D,I}) where {D,I} = D
integlevel(::BasisFunction{D,I}) where {D,I} = I

unitdegree(::BasisFunction{D,I}) where {D,I} = D+I
unitdegree(::BasisFunction{:log,I}) where {I} = I
unitdegree(::BasisFunction{:logdivT,I}) where {I} = D-1

(bf::BasisFunction{0,0})(T) = one(T)
strexp(α) = α == 1 ? "" : normal_to_super("$α")
strcoef(α) = α == 1 ? "" : "$α"
show_expr(::BasisFunction{0,0}) = ""
for α in 1:5
    @eval begin
        (bf::BasisFunction{$α,0})(T) = T^$α
        (bf::BasisFunction{-$α,0})(T) = inv(T^$α)
        show_expr(::BasisFunction{$α,0}) = "T$(strexp($α))"
        show_expr(::BasisFunction{-$α,0}) = "/T$(strexp($α))"
    end
end
(bf::BasisFunction{0.5,0})(T) = sqrt(T)
show_expr(::BasisFunction{0.5,0}) = "√T"
(bf::BasisFunction{-0.5,0})(T) = inv(sqrt(T))
show_expr(::BasisFunction{-0.5,0}) = "/√T"
(bf::BasisFunction{:log,0})(T) = log(ustrip(T))
show_expr(::BasisFunction{:log,0}) = "ln(T)"
(bf::BasisFunction{:logdivT,0})(T) = log(ustrip(T))/T
show_expr(::BasisFunction{:logdivT,0}) = "ln(T)/T"

(bf::BasisFunction{0,-1})(T) = zero(T)
for α in (1:5..., 0.5, 1.5)
    @eval begin
        (bf::BasisFunction{$α,-1})(T) = $α * T^($α - 1)
        show_expr(::BasisFunction{$α,-1}) = "$(strcoef($α))T$(strexp($α-1))"
        (bf::BasisFunction{-$α,-1})(T) = -$α * inv(T^($α + 1))
        show_expr(::BasisFunction{-$α,-1}) = "(-$(strcoef($α))/T$(strexp($α+1)))"
    end
end
(bf::BasisFunction{:log,-1})(T) = inv(T)
show_expr(::BasisFunction{:log,-1}) = "/T"
(bf::BasisFunction{:logdivT,-1})(T) = (1 - log(ustrip(T)))/T^2
show_expr(::BasisFunction{:logdivT,-1}) = "(1-ln(T))/T²"

for α in (0:5..., 0.5, 1.5)
    @eval begin
        (bf::BasisFunction{$α,1})(T) = T^($α + 1) / ($α + 1)
        if $α != 0 show_expr(::BasisFunction{$α,1}) = "T$(strexp($α+1))/$($α+1)" end
    end
    if α ∉ (0, 1)
        @eval begin
            (bf::BasisFunction{-$α,1})(T) = inv(T^($α - 1)) / (1 - $α)
            if $α != 2 show_expr(::BasisFunction{-$α,1}) = "(-1/($(strcoef($α-1))T$(strexp($α-1))))" end
        end
    end
end
show_expr(::BasisFunction{-2,1}) = "(-1/T)"
show_expr(::BasisFunction{0,1}) = "T"
(bf::BasisFunction{-1,1})(T) = log(ustrip(T))
show_expr(::BasisFunction{-1,1}) = "ln(T)"
(bf::BasisFunction{:log,1})(T) = T * (log(ustrip(T)) - 1)
show_expr(::BasisFunction{:log,1}) = "T(ln(T)-1)"
(bf::BasisFunction{:logdivT,1})(T) = (log(ustrip(T)))^2 / 2
show_expr(::BasisFunction{:logdivT,1}) = "ln(T)²/2"

for α in (0:5..., 0.5, 1.5)
    @eval begin
        (bf::BasisFunction{$α,2})(T) = T^($α + 2) / (($α + 1)*($α + 2))
        show_expr(::BasisFunction{$α,2}) = "T$(strexp($α+2))/$(($α+1)*($α+2))"
    end 
    if α ∉ (0, 1, 2)
        @eval begin
            (bf::BasisFunction{-$α,2})(T) = inv(T^($α - 2)) / ((1 - $α)*(2 - $α))
            show_expr(::BasisFunction{-$α,2}) = "(1/($(($α-1)*($α-2))T$(strexp($α-1))))"
        end 
    end
end
(bf::BasisFunction{-1,2})(T) = T*(log(ustrip(T)) - 1)
show_expr(::BasisFunction{-1,2}) = "T(ln(T)-1)"
(bf::BasisFunction{-2,2})(T) = -log(ustrip(T))
show_expr(::BasisFunction{-2,2}) = "(-ln(T))"
(bf::BasisFunction{:log,2})(T) = T^2*log(ustrip(T))/2 - 3*T^2/4
show_expr(::BasisFunction{:log,2}) = "(T²ln(T)/2 - 3/4T²)"
(bf::BasisFunction{:logdivT,2})(T) = T*(log(ustrip(T)))^2/2 - T*log(ustrip(T)) + T
show_expr(::BasisFunction{:logdivT,2}) = "(Tln(T)²/2 - Tln(T) + T)"

BasisFunction{1,-1}() = BasisFunction{0,0}()

∂(::BasisFunction{D,I}) where {D,I} = BasisFunction{D,I-1}()
∫(::BasisFunction{D,I}) where {D,I} = BasisFunction{D,I+1}()
∬(::BasisFunction{D,I}) where {D,I} = BasisFunction{D,I+2}()
divT(::BasisFunction{D,I}) where {D,I} = BasisFunction{D-1,I}()
divT(::BasisFunction{:log,I}) where {I} = BasisFunction{:logdivT,I}()

isconstant(::BasisFunction) = false
isconstant(::BasisFunction{0,0}) = true
iszerofunc(::BasisFunction) = false
iszerofunc(::BasisFunction{0,-1}) = true

show_expr(bf::BasisFunction) = "unknown"

function Base.show(io::IO, bf::BasisFunction)
    print(io, show_expr(bf))
end

struct ThermoFunction{F<:NamedTuple,C<:NamedTuple,T<:Number,Z<:Number} <: Callable
    bases::F
    coeffs::C
    Tref::T
    zeroinit::Z
    cstidx::Union{Symbol,Nothing}
    param::Symbol
end

function (lf::ThermoFunction)(T)
    s = lf.zeroinit
    @inbounds @simd for name in keys(lf.coeffs)
        s += getfield(lf.coeffs, name) * getfield(lf.bases, name)(T)
    end
    return s
end

function (lf::ThermoFunction)()
    return lf(lf.Tref)
end

coefficients(lf::ThermoFunction) = lf.coeffs

function ThermoFunction(degrees::AbstractVector, coeffs::AbstractVector{<:Number}; Tref=298.15, startindex=0, param=:a)

    if length(coeffs) == 0
        cstidx = Symbol(param, "₀")
        return ThermoFunction((cstidx => BasisFunction{0, 0}()), (cstidx => 0), Tref, 0, cstidx, param)
    end

    with_units = promote_type(typeof.(coeffs)...) <: Quantity
    if with_units && !(Tref isa Quantity)
        Tref *= K
    end

    funcs = [BasisFunction{d, 0}() for d in degrees]

    nonzero = findall(!iszero, coeffs)
    if length(nonzero) == 0
        nonzero = [1]
    end
    kept_funcs = funcs[nonzero]
    kept_coefs = coeffs[nonzero]

    varunit = unit(1)
    if with_units
        varunitvec = unit.(kept_coefs) .* (K .^ unitdegree.(kept_funcs))
        if !allequal(varunitvec) error("Degrees $degrees and coefficients $coeffs are not dimensionally compatible.") end
        varunit = length(varunitvec)>0 ? varunitvec[1] : unit(1)
    end

    kept_names = [Symbol(param, normal_to_sub(string(i + startindex - 1))) for i in nonzero]
    idxcst = findfirst(==(0), degrees)
    cstidx = isnothing(idxcst) ? nothing : Symbol(param, normal_to_sub(string(idxcst + startindex - 1)))

    base_nt = NamedTuple{Tuple(kept_names)}(Tuple(kept_funcs))
    coef_nt = NamedTuple{Tuple(kept_names)}(Tuple(kept_coefs))

    return ThermoFunction(base_nt, coef_nt, Tref, with_units ? 0*varunit : 0, cstidx, param)
end

function +(tf::ThermoFunction, x::Number)
    if isnothing(tf.cstidx)
        cstidx = Symbol(tf.param, "₀")
        idx = 0
        while hasfield(typeof(tf.coeffs), cstidx)
            idx -= 1
            cstidx = Symbol(tf.param, normal_to_sub(string(idx)))
        end
        return ThermoFunction(merge(NamedTuple{(cstidx,)}((BasisFunction{0,0}(),)), tf.bases),
                              merge(NamedTuple{(cstidx,)}(x), tf.coeffs),
                              tf.Tref, tf.zeroinit, cstidx, tf.param)
    else
        return ThermoFunction(tf.bases, merge(tf.coeffs, NamedTuple{(tf.cstidx,)}((getfield(tf.coeffs, tf.cstidx)+x,))), tf.Tref, tf.zeroinit, tf.cstidx, tf.param)
    end
end

+(x::Number, tf::ThermoFunction) = +(tf, x)

function -(tf::ThermoFunction, x::Number)
    if isnothing(tf.cstidx)
        cstidx = Symbol(tf.param, "₀")
        idx = 0
        while hasfield(typeof(tf.coeffs), cstidx)
            idx -= 1
            cstidx = Symbol(tf.param, normal_to_sub(string(idx)))
        end
        return ThermoFunction(merge(NamedTuple{(cstidx,)}((BasisFunction{0,0}(),)), tf.bases),
                              merge(NamedTuple{(cstidx,)}(-x), tf.coeffs),
                              tf.Tref, tf.zeroinit, cstidx, tf.param)
    else
        return ThermoFunction(tf.bases, merge(tf.coeffs, NamedTuple{(tf.cstidx,)}((getfield(tf.coeffs, tf.cstidx)-x,))), tf.Tref, tf.zeroinit, tf.cstidx, tf.param)
    end
end

function *(tf::ThermoFunction, x::Number)
    ThermoFunction(tf.bases, (; (k => v*x for (k,v) in pairs(tf.coeffs))...), tf.Tref, tf.zeroinit, tf.cstidx, tf.param)
end

*(x::Number, tf::ThermoFunction) = *(tf, x)

-(tf::ThermoFunction) = *(tf, -1)

-(x::Number, tf::ThermoFunction) = +(-tf, x)

function /(tf::ThermoFunction, x::Number)
    ThermoFunction(tf.bases, (; (k => v/x for (k,v) in pairs(tf.coeffs))...), tf.Tref, tf.zeroinit, tf.cstidx, tf.param)
end

function ∫(tf::ThermoFunction)
    newbases = (; (k => ∫(v) for (k,v) in pairs(tf.bases))...)
    ThermoFunction(newbases, tf.coeffs, tf.Tref, tf.Tref isa Quantity ? tf.zeroinit*K : tf.zeroinit, nothing, tf.param)
end

function ∬(tf::ThermoFunction)
    newbases = (; (k => ∬(v) for (k,v) in pairs(tf.bases))...)
    ThermoFunction(newbases, tf.coeffs, tf.Tref, tf.Tref isa Quantity ? tf.zeroinit*K^2 : tf.zeroinit, nothing, tf.param)
end

function ∂(tf::ThermoFunction)
    newbases = (; (k => ∂(v) for (k,v) in pairs(tf.bases))...)
    newcoeffs = tf.coeffs
    zeroidx = findfirst(bf->iszerofunc(bf), newbases)
    if !isnothing(zeroidx)
        newbases = (; (k => v for (k,v) in pairs(newbases) if k != zeroidx)...)
        newcoeffs = (; (k => v for (k,v) in pairs(newcoeffs) if k != zeroidx)...)
    end
    cstidx = findfirst(bf->isconstant(bf), newbases)
    ThermoFunction(newbases, newcoeffs, tf.Tref, tf.Tref isa Quantity ? tf.zeroinit/K : tf.zeroinit, cstidx, tf.param)
end

function divT(tf::ThermoFunction)
    newbases = (; (k => divT(v) for (k,v) in pairs(tf.bases))...)
    cstidx = findfirst(bf->isconstant(bf), newbases)
    ThermoFunction(newbases, tf.coeffs, tf.Tref, tf.Tref isa Quantity ? tf.zeroinit/K : tf.zeroinit, cstidx, tf.param)
end

function Base.show(io::IO, lf::ThermoFunction)
    print(io, join(["$k$v" for (k,v) in pairs(lf.bases)]," + "))
    print(io, " with {")
    print(io, replace(string(lf.coeffs), "("=>"", ")"=>"", " = "=>"="))
    print(io, " ; Tref=", lf.Tref,"}")
end
