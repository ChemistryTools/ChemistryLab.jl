using Symbolics
using DynamicQuantities

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

# function extract_vars(expr)
#     vars = Set{Symbol}()
#     function _extract(e)
#         if isa(e, Symbol)
#             push!(vars, e)
#         elseif isa(e, Expr)
#             for arg in e.args
#                 _extract(arg)
#             end
#         end
#     end
#     _extract(expr)
#     return collect(vars)
# end

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

Cp = :(a₀ + a₁*T + a₂/T^2 + a₃/√T + a₄*T^2 + a₅*T^3 + a₆*T^4 + a₇/T^3 + a₈/T + a₉*√T + a₁₀*log(T))

v = sort(extract_vars(Cp); by = x->Int(parse_symbol_num(string(x)).index))
deleteat!(v, findall(x -> x == :T, v))

dict_vars = Dict{Symbol, Num}()
for x in v @eval dict_vars[Symbol($x)] = (@variables $x)[1] end

degrees = [0, 1, -2, -0.5, 2, 3, 4, -3, -1, 0.5, :log]
coeffs = [210.0J/K/mol, 0.120J/mol/K^2, -3.07e6J*K/mol, 0.0J/mol/√K]

params = NamedTuple{Tuple(v[CartesianIndices(coeffs)])}(Tuple(coeffs))
params_nounit = NamedTuple{keys(params)}(ustrip.(values(params)))
subs = Dict(p for p in pairs(params_nounit))
for (k,v) in dict_vars
    if !haskey(subs, k) && k != :T
        subs[k] = 0
    end
end
# numCp2(T) = eval(substitute_expr(Cp, subs))

numCp = eval(substitute_expr(Cp, dict_vars))

Cpr = substitute(eval(Cp), Dict(eval(Symbol(:a,subscriptnumber(i)))=>0 for i in 4:10))

Cpf = eval(build_function(substitute(Cpr, Dict(eval(k)=>v for (k,v) in subs)), T))


# struct ThermoFunction{N,V,P,U,F,FU}
#     symbolic::N
#     variables::V
#     parameters::P
#     unit::U
#     func::F
#     func_unit::FU
# end

# function ThermoFunction(expr)

#     @variables T, a₀, a₁, a₂
#     symbolic = a₀ + a₁ * T + a₂ / T^2
#     func = eval(build_function(symbolic, T, a₀, a₁, a₂; expression=Val(true)))
#     return ThermoFunction(symbolic, func)
# end

# function deriv(expr, var)
#     return Symbolics.derivative(expr, (@eval @variables var)[1])
# end

# function derivT(expr)
#     @variables T
#     return Symbolics.derivative(expr, T)
# end

# function (tf::ThermoFunction)(args...)
#     return tf.func(args...)
# end



