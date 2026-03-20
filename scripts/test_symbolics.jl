using DynamicQuantities
using ModelingToolkit
using SymbolicNumericIntegration

function subscriptnumber(i::Integer)
    if i < 0
        c = [Char(0x208B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        push!(c, Char(0x2080+d))
    end
    return join(c)
end

function superscriptnumber(i::Integer)
    if i < 0
        c = [Char(0x207B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end

function from_subscriptnumber(s::String)
    chars = collect(s)
    negative = !isempty(chars) && chars[1] == Char(0x208B)
    if negative
        chars = chars[2:end]
    end
    value = 0
    for c in chars
        digit = Int(c) - 0x2080
        value = value * 10 + digit
    end
    return negative ? -value : value
end

function from_superscriptnumber(s::String)
    chars = collect(s)
    negative = !isempty(chars) && chars[1] == Char(0x207B)
    if negative
        chars = chars[2:end]
    end
    value = 0
    for c in chars
        if c == Char(0x2070)
            digit = 0
        elseif c == Char(0x00B9)
            digit = 1
        elseif c == Char(0x00B2)
            digit = 2
        elseif c == Char(0x00B3)
            digit = 3
        else
            digit = Int(c) - 0x2070
        end
        value = value * 10 + digit
    end
    return negative ? -value : value
end

function parse_symbol_num(s::AbstractString)
    chars = collect(s)
    n = length(chars)
    pos = n + 1
    for i in 1:n
        c = chars[i]
        if '0' <= c <= '9'
            pos = i
            break
        end
        if 0x2080 <= Int(c) <= 0x2089
            pos = i
            break
        end
        if c in (Char(0x2070), Char(0x00B9), Char(0x00B2), Char(0x00B3)) || (0x2074 <= Int(c) <= 0x2079)
            pos = i
            break
        end
    end
    if pos == n + 1
        return (base=s, index=-100, convert_func=identity)
    end
    radical = join(chars[1:pos-1])
    number_chars = chars[pos:end]
    firstnum = number_chars[1]
    number = join(number_chars)
    if '0' <= firstnum <= '9'
        nature = identity
        number = parse(Int,number)
    elseif 0x2080 <= Int(firstnum) <= 0x2089
        nature = subscriptnumber
        number = from_subscriptnumber(number)
    elseif firstnum in (Char(0x2070), Char(0x00B9), Char(0x00B2), Char(0x00B3)) || (0x2074 <= Int(firstnum) <= 0x2079)
        nature = superscriptnumber
        number = from_superscriptnumber(number)
    else
        nature = identity
        number = 0
    end
    return (base=radical, index=number, convert_func=nature)
end

function extract_vars(expr)
    vars = Set{Symbol}()
    function _extract(e, is_func=false)
        if isa(e, Symbol)
            if !is_func
                push!(vars, e)
            end
        elseif isa(e, Expr)
            if e.head == :call
                _extract(e.args[1], true)
                for arg in e.args[2:end]
                    _extract(arg, false)
                end
            else
                for arg in e.args
                    _extract(arg, false)
                end
            end
        end
    end
    _extract(expr, false)
    return collect(vars)
end

function substitute_expr(expr, subs::AbstractDict)
    for (k, v) in subs
        if isequal(expr, k)
            return v
        end
    end
    if expr isa Expr
        return Expr(expr.head, (substitute_expr.(expr.args, Ref(subs)))...)
    else
        return expr
    end
end

function symbolic_to_expr(exp)
    if exp isa Num && isnumeric(exp)
        return :( $(exp) )
    elseif exp isa Num
        args_expr = map(symbolic_to_expr, arguments(exp))
        return Expr(:call, operation(exp), args_expr...)
    else
        return exp
    end
end

struct SymbolicFunc{P,U,F,R}
    expr::Num
    parameters::P
    var::Num
    unit::U
    func::F
    ref::R
end

function SymbolicFunc(symexpr, params_unit::NamedTuple, var=nothing; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    symexpr = Num(symexpr)
    varofexpr = Symbol.(get_variables(symexpr))
    if isnothing(var)
        vars = [:T, :t]
        idxvar = findfirst(x->x in varofexpr, vars)
        var = @eval (@variables $(vars[idxvar]))[1]
    end
    # existvars = in.(collect(keys(params_unit)), Ref(varofexpr))
    # params_unit = NamedTuple{keys(params_unit)[existvars]}(values(params_unit)[existvars])
    params_unit = (; (k => v for (k,v) in pairs(params_unit) if k in varofexpr)...)
    # params_nounit = NamedTuple{keys(params_unit)}(ustrip.(values(params_unit)))
    params_nounit = (; (k=>ustrip(v) for (k,v) in pairs(params_unit))...)
    func = eval(build_function(substitute(symexpr, [eval(k)=>v for (k,v) in pairs(params_nounit)]), var))
    nounitfunc = [f(var) => 1 for f in [log, log10, log2, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, exp]]
    unit = promote_type(typeof.(values(params_unit))...) <: Quantity ?
                   Quantity(1, dimension(substitute(symexpr, [[eval(k)=>v for (k,v) in pairs(params_unit)] ; nounitfunc ; var=>getfield(ref, Symbol(var))]))) : 1
    return SymbolicFunc(symexpr, params_unit, var, unit, func, ref)
end

function SymbolicFunc(expr::Expr, params::AbstractVector{<:Number}, var=nothing; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    v = extract_vars(expr)
    if isnothing(var)
        vars = [:T, :t]
        idxvar = findfirst(x->x in v, vars)
        var = @eval (@variables $(vars[idxvar]))[1]
    end
    deleteat!(v, findfirst(x -> x == Symbol(var), v))
    sort!(v; by = x->Int(parse_symbol_num(string(x)).index))
    dict_vars = Dict{Symbol, Num}()
    for x in v dict_vars[x] = @eval (@parameters $x)[1] end
    subszeros = Dict(eval(k)=>0 for (i,k) in enumerate(v) if i>length(params) || iszero(params[i]))
    symexpr = substitute(eval(expr), subszeros)
    nonzeros = (!iszero).(params)
    params_unit = NamedTuple{Tuple(v[CartesianIndices(params)][nonzeros])}(Tuple(params[nonzeros]))
    return SymbolicFunc(symexpr, params_unit, var; ref=ref)
end

const thermo_function_library = Dict(
    :Cp => :(a₀ + a₁*T + a₂/T^2 + a₃/√T + a₄*T^2 + a₅*T^3 + a₆*T^4 + a₇/T^3 + a₈/T + a₉*√T + a₁₀*log(T)),
    :logKr => :(A₀ + A₁*T + A₂/T + A₃*log(T) + A₄/T^2 + A₅*T^2 + A₆*√T)
)

function SymbolicFunc(sym::Symbol, params::AbstractVector{<:Number}, var=nothing; ref=(T=298.15u"K", P=1u"bar", t=0u"s"))
    return SymbolicFunc(thermo_function_library[sym], params, var; ref=ref)
end

function (tf::SymbolicFunc)(T)
    return tf.func(T)
end

function (tf::SymbolicFunc)(T::Quantity)
    return tf.func(ustrip(T))*tf.unit
end

import Base: +, -, *, /

function +(tf::SymbolicFunc, x::Number)
    indexcst = findfirst(a->iszero(expand_derivatives(Differential(tf.var)(Differential(eval(a))(tf.expr)))), keys(tf.parameters))
    if isnothing(indexcst)
        cst = :cst
        if length(tf.parameters) > 0
            ps = parse_symbol_num(string(keys(tf.parameters)[1]))
            cst = Symbol(ps.base, ps.convert_func(ps.index-1))
        end
        @eval @parameters $(cst)
        expr = eval(cst) + tf.expr
        return SymbolicFunc(expand(expr), merge((NamedTuple{(cst,)}(x)), tf.parameters), tf.var; ref=tf.ref)
    else
        cst = keys(tf.parameters)[indexcst]
        return SymbolicFunc(tf.expr, merge(tf.parameters, NamedTuple{(cst,)}((getfield(tf.parameters, cst)+x,))), tf.var; ref=tf.ref)
    end
end

+(x::Number, tf::SymbolicFunc) = +(tf, x)

-(tf::SymbolicFunc, x::Number) = +(tf, -1*x)

-(tf::SymbolicFunc) = SymbolicFunc(-1*tf.expr, tf.parameters, tf.var; ref=tf.ref)

*(tf::SymbolicFunc, x::Number) = SymbolicFunc(x*tf.expr, tf.parameters, tf.var; ref=tf.ref)

*(x::Number, tf::SymbolicFunc) = *(tf, x)

/(tf::SymbolicFunc, x::Number) = SymbolicFunc(tf.expr/x, tf.parameters, tf.var; ref=tf.ref)

∂(tf::SymbolicFunc, var=tf.var) = SymbolicFunc(expand(expand_derivatives(Differential(var)(tf.expr))), tf.parameters, tf.var; ref=tf.ref)

∫(tf::SymbolicFunc, var=tf.var) = SymbolicFunc(expand(expand_derivatives(integrate(tf.expr, var; symbolic=true, detailed=false))), tf.parameters, tf.var; ref=tf.ref)

function Base.show(io::IO, tf::SymbolicFunc)
    print(io, tf.expr)
    print(io, " with {")
    print(io, replace(string(tf.parameters), "("=>"", ")"=>"", " = "=>"="))
    print(io, " ; ref: ", replace(string(tf.ref), "("=>"", ")"=>"", " = "=>"="), "}")
end

expr = :(a₀ + a₁*T + a₂/T^2 + a₃/√T + a₄*T^2 + a₅*T^3 + a₆*T^4 + a₇/T^3 + a₈/T + a₉*√T + a₁₀*log(T))
params = [210.0u"J/K/mol", 0.0u"J/mol/K^2", -3.07e6u"J*K/mol", 0.10u"J/mol/√K"]

tf = SymbolicFunc(expr, params)

tf2 = SymbolicFunc(expand(expand_derivatives(Differential(T)(tf.expr))), tf.parameters)

tf = SymbolicFunc(:((c₁+c₂*t)/(c₃+c₄*t^2)), [1u"J",2u"J/s",3,4u"1/s^2"]; ref=(t=0u"s",))

Cp = SymbolicFunc(dict_cp_ft_equation[:Cp], params)
