@testset "Simplex" begin
    @testset "Primal" begin
        @testset "Two Variable" begin
            A = Float64[-1 1 1 0 0 11
                1 1 0 1 0 27
                2 5 0 0 1 90
                -4 -6 0 0 0 0]
            answer = simplex(A)[2]
            true_answer = Dict{Any,Any}("x_1" => 15, "x_2" => 12)
            for (key, value) in true_answer
                @test answer[key] ≈ value
            end
        end

        @testset "Three Variable" begin
            A = Float64[1 2 3/2 1 0 0 12
                2/3 2/3 1 0 1 0 4.6
                1/2 1/3 1/2 0 0 1 2.4
                -11 -16 -15 0 0 0 0]
            B = Float64[4 1 1 1 0 0 30
                2 3 1 0 1 0 60
                1 2 3 0 0 1 40
                -3 -2 -1 0 0 0 0]
            data = Dict{Any,Any}("A" => A, "B" => B)
            true_answer = Dict{Any,Any}(
                "A" => Dict{Any,Any}("x_1" => 0.6, "x_2" => 5.1, "x_3" => 0.8),
                "B" => Dict{Any,Any}("x_1" => 3, "x_2" => 18),
            )
            for (name, array) in data
                ans = simplex(array)[2]
                for (key, answer) in true_answer[name]
                    @test ans[key] ≈ answer
                end
            end
        end
    end
    @testset "Dual" begin
        A = Float64[60 12 10 1 0 0.12
            60 6 30 0 1 0.15
            -300 -36 -90 0 0 0]
        true_answer = Dict{Any,Any}("x_1" => 3, "x_2" => 2)
        ans = simplex(A; primal=false)[2]
        for (key, answer) in true_answer
            @test ans[key] ≈ answer
        end

    end

end

@testset "simplex_case" begin
    @testset "min_base" begin
        A = Float64[60 60
            12 6
            10 30]
        b = Float64[300; 36; 90]
        c = Float64[0.12; 0.15; 0]
        inequality = [">=" for i in 1:3]
        answer = simplex_case(A, b, c; type="min_base", inequality=inequality)[2]
        true_answer = Dict{Any,Any}("x_1" => 3, "x_2" => 2)
        for (key, value) in true_answer
            @test answer[key] ≈ value
        end
    end
    @testset "max_base" begin
        A = Float64[-1 1
            1 1
            2 5]
        b = Float64[11; 27; 90]
        c = Float64[4; 6; 0]
        answer = simplex_case(A, b, c; type="max_base")[2]
        true_answer = Dict{Any,Any}("x_1" => 15, "x_2" => 12)
        for (key, value) in true_answer
            @test answer[key] ≈ value
        end
    end
    @testset "max_mixed" begin
        A = Float64[10 5
            2 3
            1 0
            0 1]
        b = Float64[200; 60; 12; 6]
        c = Float64[2000; 1000; 0]
        inequality = ["<=", "=", "<=", ">="]
        answer = simplex_case(A, b, c; type="max_mixed", inequality=inequality)[2]
        true_answer = Dict{Any,Any}("x_1" => 12, "x_2" => 12)
        for (key, value) in true_answer
            @test answer[key] ≈ value
        end
    end
end
"Done"
