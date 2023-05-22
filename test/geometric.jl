
const R = OpenLocationCode.GEO_RADIUS

@testset "distances" begin
    lat1, lon1 = 50.1732, 8.3381
    lat2, lon2 = 50.0279, 8.0334
    c1 = encode(lat1, lon1, 15)
    c2 = encode(lat2, lon2, 15)
    @test distance(c1, c2) ≈ 27111.7984
    @test distance(-90.0, 350, 0, 0) ≈ π / 2 * R
end

@testset "areas" begin
    a00 = area(encode(-10, 0, 2))
    a10 = area(encode(10, 0, 2))
    a30 = area(encode(30, 0, 2))
    a50 = area(encode(50, 0, 2))
    a70 = area(encode(70, 0, 2))
    s20 = a00 + 2*(a10 + a30 + a50 + a70)
    s00 = s20 * 18
    @test s00 ≈ 4π *  R ^ 2
end
