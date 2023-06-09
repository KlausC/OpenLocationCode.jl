# Test data for validity tests.
# Format of each line is:
#   code,isValid,isShort,isFull
# Valid full codes:
DATA = [
("8FWC2345+G6",true,false,true)
("8FWC2345+G6G",true,false,true)
("8fwc2345+",true,false,true)
("8FWCX400+",true,false,true)
# Valid short codes:
("WC2345+G6g",true,true,false)
("2345+G6",true,true,false)
("45+G6",true,true,false)
("+G6",true,true,false)
# Invalid codes
("G+",false,false,false)
("+",false,false,false)
("8FWC2345+G",false,false,false)
("8FWC2_45+G6",false,false,false)
("8FWC2η45+G6",false,false,false)
("8FWC2345+G6+",false,false,false)
("8FWC2345G6+",false,false,false)
("8FWC2300+G6",false,false,false)
("WC2300+G6g",false,false,false)
("WC2345+G",false,false,false)
("WC2300+",false,false,false)
# Validate that codes at and exceeding 15 digits are still valid when all their
# digits are valid, and invalid when not.
("849VGJQF+VX7QR3J",true,false,true)
("849VGJQF+VX7QR3U",false,false,false)
("849VGJQF+VX7QR3JW",true,false,true)
("849VGJQF+VX7QR3JU",false,false,false)
("62H00002+",false,false,false)
("620AA000+",false,false,false)
("+",false,false,false)
]

@testset "validity $(d[1])" for d in DATA
    code, isValid, isShort, isFull = d
    @test is_valid(code) == isValid
    @test is_short(code) == isShort
    @test is_full(code) == isFull
end
