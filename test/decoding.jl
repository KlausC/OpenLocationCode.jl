# Test decoding Open Location Codes.
#
# Provides test cases for decoding valid codes.
#
# Format:
#   code,length,latLo,lngLo,latHi,lngHi
DATA = [
("7FG49Q00+",6,20.35,2.75,20.4,2.8)
("7FG49QCJ+2V",10,20.37,2.782125,20.370125,2.78225)
("7FG49QCJ+2VX",11,20.3701,2.78221875,20.370125,2.78225)
("7FG49QCJ+2VXGJ",13,20.370113,2.782234375,20.370114,2.78223632813)
("8FVC2222+22",10,47.0,8.0,47.000125,8.000125)
("4VCPPQGP+Q9",10,-41.273125,174.785875,-41.273,174.786)
("62G20000+",4,0.0,-180.0,1,-179)
("22220000+",4,-90,-180,-89,-179)
("7FG40000+",4,20.0,2.0,21.0,3.0)
("22222222+22",10,-90.0,-180.0,-89.999875,-179.999875)
("6VGX0000+",4,0,179,1,180)
("6FH32222+222",11,1,1,1.000025,1.00003125)
("CFX30000+",4,89,1,90,2)
("62H20000+",4,1,-180,2,-179)
("62H30000+",4,1,-179,2,-178)
("CFX3X2X2+X2",10,89.9998750,1,90,1.0001250)
# Test non-precise latitude/longitude value
("6FH56C22+22",10,1.2000000000000028,3.4000000000000057,1.2001249999999999,3.4001250000000027)
# Validate that digits after the first 15 are ignored when decoding
("849VGJQF+VX7QR3J",15,37.5396691200,-122.3750698242,37.5396691600,-122.3750697021)
("849VGJQF+VX7QR3J7QR3J",15,37.5396691200,-122.3750698242,37.5396691600,-122.3750697021)
]

@testset "decoding $(d[1])" for d in DATA
    code, codelength, latlo, longlo, lathi, longhi = d
    ca = decode(code)
    @test ca ≈ CodeArea(latlo, longlo, codelength)
    @test latitude_high(ca) ≈ lathi
    @test longitude_high(ca) ≈ longhi
end

@testset "validity" begin
    @test  is_full("C2XXXXXX+")
    @test !is_full("F2XXXXXX+")
    @test  is_full("2VXXXX00+")
    @test !is_full("2WXXXXXX+")
    @test !is_valid("+0022")
    @test !is_valid("00+")
    @test_throws ArgumentError decode("XX+22")
    @test_throws ArgumentError recover_nearest("X+", 0, 0)
    @test_throws ArgumentError shorten("XX+22", 0, 0)
    @test_throws ArgumentError shorten("22XX0000+", 0, 0)
end

@testset "show and convert CodeArea" begin
    code = "9G000000+"
    ca = decode(code)
    @test latlong(ca) == (60.0, 30.0)
    @test sprint(show, ca) == "CodeArea(\"$code\", 50.0+20.0, 20.0+20.0, 2)"
    @test CodeArea(code) == ca
    @test CodeArea(50, 20, 2) == ca
end
