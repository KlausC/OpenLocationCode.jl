# OpenLocationCode

[![Build Status](https://github.com/KlausC/OpenLocationCode.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KlausC/OpenLocationCode.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Convert locations to and from short codes

Open Location Codes or "Plus Codes" are short, 10-15 character codes that can be used instead
of street addresses. The codes can be generated and decoded offline, and use
a reduced character set that minimises the chance of codes including words.

Codes are able to be shortened relative to a nearby location. This means that
in many cases, only four to seven characters of the code are needed.
To recover the original code, the same location is not required, as long as
a nearby location is provided.

Codes represent rectangular areas rather than points, and the longer the
code, the smaller the area. A 10 character code represents a 13.9x13.9
meter area (at the equator). An 11 character code represents approximately
a 2.8x3.5 meter area, while a 15 character code is 4.5x13.6 mm.

Two encoding algorithms are used. The first 10 characters are pairs of
characters, one for latitude and one for longitude, using base 20. Each pair
reduces the area of the code by a factor of 400. Only even code lengths are
sensible, since an odd-numbered length would have sides in a ratio of 20:1.

At positions 11-15, the algorithm changes so that each character selects one
position from a 4x5 grid. This allows refinements by one to five characters.

Examples:

```julia
julia> using OpenLocationCode

julia> # Encode a location, default accuracy:
       encode(47.365590, 8.524997)
"8FVC9G8F+6X"

julia> # Encode a location using one stage of additional refinement:
       encode(47.365590, 8.524997, 15)
"8FVC9G8F+6XQQ435"

julia> # Decode a full code:
       code = "8FVCCJ8F+6X"
"8FVCCJ8F+6X"

julia> ca = decode(code)
CodeArea{Float64}(47.4155, 8.624875, 10)

julia> println("# Center is lat=$(latitude_center(ca)), lon=$(longitude_center(ca))")
# Center is lat=47.4155625, lon=8.6249375

julia> println("# extension of area is $(latitude_precision(ca) * 111321)x$(longitude_precision(ca) * 111321 * cosd(latitude_low(ca))) m")
# extension of area is 13.915125x9.416042416499675 m

julia> # Attempt to trim the first characters from a code:
       shorten("8FVC9G8F+6X", 47.5, 8.5)
"9G8F+6X"

julia> # Recover the full code from a short code:
       recover_nearest("9G8F+6X", 47.4, 8.6)
"8FVC9G8F+6X"

julia> recover_nearest("8F+6X", 47.4, 8.6)
"8FVCCJ8F+6X"
```

This `Julia` implementation originates from the `Python` sources and test data
in [google/open-location-code](https://github.com/google/open-location-code) v1.0.4
with this [specification](https://github.com/google/open-location-code/blob/main/docs/specification.md).

© 2023 Klaus Crusius
