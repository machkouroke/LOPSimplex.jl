include("utilities.jl")
"""
    detect_solution(simplex_array::Matrix{Float64}, solution_set::Vector{String}; primal=true)
Permet d'extraire les solutions aux problème. Pour le problème primal les solutions sont des variables de bases tandis que pour le problèmes
dual les solutions sont des variables d'écart
# Arguments:
- `simplex_array::Matrix{Float64}`: Matrice final du tableau de simplex	
- `solution_set::Vector{String}`: Pour le problème primal, les variable de bases sont envoyé. 
Pour le problème dual, toutes les variable sont envoyé
- `primal::Bool`: Si `true` alors le problème est primal, sinon le problème est dual
"""
function detect_solution(simplex_array::Matrix{Float64}, solution_set::Vector{String}; primal=true)
    variables = Dict()
    if primal
        for (index, variable) in enumerate(solution_set)
            if split(variable, "_")[1] == "x"
                variables[variable] = simplex_array[index, end]
            end
        end
    else
        for (index, variable) in enumerate(solution_set)
            if split(variable, "_")[1] == "e"
                variables["x_$(split(variable, "_")[2])"] = simplex_array[end, index]
            end
        end
    end
    return variables
end

"""
    incoming(x::Matrix{Float64})
Permet de trouver la variable entrante dans la base. Vu qu'un signe - est affecté aux coefficients de la fonction objectif, celle ci seras celle
qui a le plus petit coefficient négatif.
# Arguments:
- `x::Matrix{Float64}`: Matrice du tableau de simplex
"""
function incoming(x::Matrix{Float64})::Int64
    return findmin(x[end, 1:end-1])[2]
end


"""
    incoming(x::Matrix{Float64})
Permet de trouver la variable sortante de la base en utilisant la règle de Bland. 
# Arguments:
- `x::Matrix{Float64}`: Matrice du tableau de simplex
- `pivot::Int64`: Indice de la variable entrante
"""
function outgoing(x::Matrix{Float64}, pivot::Int64)::Int64
    x_div_pivot = x[1:end-1, end] ./ x[1:end-1, pivot]
    x_div_pivot_positif = ifelse.(x_div_pivot .>= 0, x_div_pivot, Inf)
    return findmin(x_div_pivot_positif)[2]
end

"""
    simplex(A::Matrix{Float64})
Effectue les opérations de pivot de l'algorithme du simplex sur la matrice A.
# Arguments:
- `A::Matrix{Float64}`: Matrice du simplex à pivoter.
- `in_base::Vector{String}`: Variable de bases dans la matrice de départ du simplex
- `all_base::Vector{String}`:Toutes les variable du système
- `verbose::Bool`: Précise si l'on souhaite afficher les étapes de l'algorithme du simplex.
- `primal::Bool`: Précise si l'on souhaite afficher la solution du problème primal.
"""
function simplex(A::Matrix{Float64}; in_base=Nothing, all_base=Nothing, verbose::Bool=false, primal=true)
    n::Int64, m::Int64 = size(A)[1] - 1, size(A)[2] - 1
    if in_base == Nothing
        in_base = ["e_$(i)" for i in 1:n]
    end
    if all_base == Nothing
        out_base::Vector{String} = ["x_$(i)" for i in 1:m-n]
        all_base = vcat(out_base, in_base)
    end
    B::Matrix{Float64} = deepcopy(A)

    while any(i -> i < 0, B[end, 1:end-1])
        verbose && @show in_base
        verbose && display(B)
        k = incoming(B)
        p = outgoing(B, k)
        in_base[p] = all_base[k]

        not_outgoing = setdiff(1:size(B)[1], [p])
        B[p, :] = B[p, :] ./ B[p, k]
        B[not_outgoing, :] = B[not_outgoing, :] - B[not_outgoing, k] * B[p, :]'
    end
    verbose && println("Final solution")
    verbose && @show in_base
    verbose && display(B)
    return detect_solution(B, primal ? in_base : all_base, primal=primal)
end



"""
    simplex_case(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}; inequality=["<=" for i in 1:size(A)[1]], type="max_base")
Permet de détecter quelle algorithme du simplex appliquer en fonction du système proposé
# Arguments:
- `A::Matrix{Float64}`: Matrice des coefficients du système de contraintes
- `b::Vector{Float64}`: Vecteur des valeurs des contraintes
- `c::Vector{Float64}`: Vecteur des coefficients de la fonction objectif
- `inequality::Vector{String}`: Vecteur des inégalités du système de contraintes
- `type::String`: Type de problème à résoudre soit une maximisation ou une minimisation
"""
function simplex_case(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64};
    inequality=["<=" for i in 1:size(A)[1]], type="max_base")
    @match type begin
        "max_base" => begin
            simp_array = simplex_matrix_builder(A, b, -c)
            answer = simplex(simp_array[1]; in_base=simp_array[3], all_base=simp_array[2], verbose=true)
        end
        "min_base" => begin
            simp_array = simplex_matrix_builder(convert(Matrix{Float64}, A'), c[1:end-1], [-b; 0])
            answer = simplex(simp_array[1]; in_base=simp_array[3], all_base=simp_array[2], verbose=true, primal=false)
        end
    end
end

function main()
    A = Float64[1 2
        3 2
        3 1]
    b = Float64[60; 120; 90]
    c = Float64[100; 100; 0]
    println("***Start***")
    answer = simplex_case(A, b, c; type="min_base")
    print(answer)
end

main()
