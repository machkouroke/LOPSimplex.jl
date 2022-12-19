using Match

"""
    add_slack_variable(A::Matrix{Float64}; inequality=["<=" for i in 1:size(A)[1]])
Ajout les variables d'écart à la matrix `A` du programme
# Arguments:
- `A::Matrix{Float64}`: Matrice A du problème
- `inequality`: Vecteur des inégalités du système de contraintes
"""
function add_slack_variable(A::Matrix{Float64}; inequality=["<=" for i in 1:size(A)[1]])
    n::Int64, m::Int64 = size(A)[1], size(A)[2]
    slack_variable::Matrix{Float64} = zeros(n, n)
    for (index, inequality) in enumerate(inequality)
        slack_variable[index, index] = @match inequality begin
            "<=" => 1
            ">=" => -1
            "=" => 0
        end
    end
    return hcat(A, slack_variable)
end


"""
    add_artificial_variable(A::Matrix{Float64}; inequality=["<=" for i in 1:size(A)[1]])
Ajout les variables artificielle à la matrix `A` du programme
# Arguments:
- `A::Matrix{Float64}`: Matrice A du problème
- `inequality`: Vecteur des inégalités du système de contraintes
"""
function add_artificial_variable(A::Matrix{Float64}; inequality=["<=" for i in 1:size(A)[1]])
    n::Int64, m::Int64 = size(A)[1], size(A)[2]
    slack_variable::Matrix{Float64} = zeros(n, n)
    for (index, inequality) in enumerate(inequality)
        slack_variable[index, index] = inequality in ["=", ">="]
    end
    return hcat(A, slack_variable)
end


"""
    variable_name_builder(A::Matrix{Float64})
Crée les noms de toutes les variables du système
# Arguments:
- `A::Matrix{Float64}`: Matrice A du problème
"""
function variable_name_builder(A::Matrix{Float64})
    variable_name::Vector{String} = ["x_$i" for i in 1:size(A)[2]]
    slack_variable_name::Vector{String} = ["e_$i" for i in 1:size(A)[1]]
    artificial_variable::Vector{String} = ["a_$i" for i in 1:size(A)[1]]
    return [variable_name; slack_variable_name; artificial_variable]
end


"""
    in_base_finder(inequality)::Vector{String}
Renvoie les noms des variable de base de départ du système
# Arguments:
- `A::Matrix{Float64}`: Inéquation du système de contraintes
"""
function in_base_finder(inequality)::Vector{String}
    in_base::Vector{String} = []
    for (index, inequality) in enumerate(inequality)
        variable = @match inequality begin
            "<=" => "e_$index"
            ">=" => "a_$index"
            "=" => "a_$index"
        end
        push!(in_base, variable)
    end
    return in_base
end


"""
    simplex_matrix_builder(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}; inequality=["<=" for i in 1:size(A)[1]])
Construit le tableau initial du simplex
# Arguments:
- `A::Matrix{Float64}`: Inéquation du système de contraintes
- `b::Vector{Float64}`: Vecteur des valeurs des contraintes
- `c::Vector{Float64}`: Vecteur des coefficients de la fonction objectif
- `inequality`: Vecteur des inégalités du système de contraintes
"""
function simplex_matrix_builder(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}; inequality=["<=" for i in 1:size(A)[1]])
    n::Int64, m::Int64 = size(A)[1], size(A)[2]
    variable_name = variable_name_builder(A)
    in_base = in_base_finder(inequality)
    A = A |>
        x -> add_slack_variable(x; inequality=inequality) |>
             x -> add_artificial_variable(x; inequality=inequality)


    c = vcat(c[1:end-1], zeros(n * 2), c[end])

    simplex_array = vcat(hcat(A, b), c')
    not_null_columns = vec(mapslices(col -> any(col .!= 0), simplex_array, dims=1))

    simplex_array = simplex_array[:, not_null_columns]
    variable_name = variable_name[not_null_columns[1:end-1]]
    return simplex_array, variable_name, in_base
end


"""
    function_by_artificial(A::Matrix{Float64}, in_base, all_variable)
Remplace la fonction objectif d'un tableau du simplex par la fonction somme des variables artificielles. 
Ces variable seront ensuite exprimé par les variables hors bases
# Arguments
- `A::Matrix{Float64}`: tableau du simplex de la première phase
- `in_base::Vector{String}`: vecteur des variables dans la base
- `all_variable::Vector{String}`: vecteur de toutes les variables
"""
function function_by_artificial(A::Matrix{Float64}, in_base, all_variable)
    A_copy = copy(A)
    artificial_row = findall(x -> x[1:2] == "a_", in_base)
    A_copy[end, :] = -sum(A[artificial_row, :], dims=1)
    in_base_indice = findall(x -> x in in_base, all_variable)
    A_copy[end, in_base_indice] .= 0
    return A_copy
end




function remove_artificial_column(A::Matrix{Float64}, variable_name::Vector{String})
    artificial_variable = findall(x -> x[1:2] == "a_", variable_name)
    A = A[:, setdiff(1:size(A)[2], artificial_variable)]
    variable_name = variable_name[setdiff(1:size(A)[2], artificial_variable)]
    return A, variable_name
end
