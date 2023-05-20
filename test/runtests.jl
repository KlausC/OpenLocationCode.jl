using OpenLocationCode
using Test

@testset "OpenLocationCode.jl" begin
    @testset "validity" begin include("validityTests.jl") end
    @testset "decode" begin include("decoding.jl") end
    @testset "encode" begin include("encoding.jl") end
    @testset "shorten" begin include("shortCodeTests.jl") end
end
