# OpenLocationCode

A package for Open Location Codes aka Plus Codes to represent geographical coordinates

[![Build Status][gha-img]][gha-url]    [![Coverage Status][codecov-img]][codecov-url]

Encode and decode, shorten and recover relative to a position.

This `Julia` implementation originates from the `Python` sources and test data
in [google/open-location-code](https://github.com/google/open-location-code) v1.0.4
with this [specification](https://github.com/google/open-location-code/blob/main/docs/specification.md).
A concise description is found here [OLC](https://en.wikipedia.org/wiki/Open_Location_Code).

## Description

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

All degrees are integer or floating point numbers. The created output degrees are `Float64`.
The latitude and longitude should be [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84) values.

## API

```doc
       encode(latitude, longitude[, codelength])

  Encode a location into an Open Location Code. Produces a code of the specified length, or the default length
  if no length is provided.

  The length determines the accuracy of the code. The default length is 10 characters, returning a code of
  approximately 13.9x13.9 meters. Longer codes represent smaller areas, but lengths > 15 are sub-centimetre and
  so 11 or 12 are probably the limit of useful codes.

  Arguments:
    •  latitude: A latitude in signed degrees. Will be clipped to the range -90 to 90.
    •  longitude: A longitude in signed degrees. Will be normalised to the range -180 to 180.
    •  codelength: The number of significant digits in the output code, not including any separator
       characters.

  Examples:

  julia> encode(50.173168, 8.338086, 11)
  "9F2C58FQ+768"
```

```doc

       decode(code)

  Decode an Open Location Code into the location coordinates.

  Arguments:
    •  code: The Open Location Code to decode.

  Returns:
  A CodeArea object that provides the latitude and longitude of two of the corners of the area, the center, and the length of the original code.
```

```doc
       shorten(code, latitude, longitude)

  Remove characters from the start of an OLC code. This uses a reference location to determine how many initial characters can be removed from the OLC code. The
  number of characters that can be removed depends on the distance between the code center and the reference location.

  The minimum number of characters that will be removed is four. If more than four characters can be removed, the additional characters will be replaced with the
  padding character. At most eight characters will be removed. The reference location must be within 50% of the maximum range. This ensures that the shortened code
  will be able to be recovered using slightly different locations.

  Arguments
    •  code: A full, valid code to shorten.
    •  latitude: A latitude, in signed degrees, to use as the reference point.
    •  longitude: A longitude, in signed degrees, to use as the reference point.

  Returns:
  Either the original code, if the reference location was not close enough, or the shortest code which can be used to recover from the reference location.
```

```doc
       recover_nearest(code, latitude, longitude)

  Recover the nearest matching code to a specified location. Given a short code of between four and seven characters, this recovers the nearest matching full code
  to the specified location.

  Arguments:
    •  code: A valid OLC character sequence.
    •  latitude: The latitude (in signed degrees) to use to find the nearest matching full code.
    •  longitude`: The longitude (in signed degrees) to use to find the nearest matching full code.

  Returns:
  The nearest full Open Location Code to the reference location that matches the short code. If the passed code was not a valid short code, but was a valid full code, it is returned with proper capitalization but otherwise unchanged.
```

## Examples

```jldoctest
julia> using OpenLocationCode

julia> # Encode a location, default accuracy:
       encode(47.365590, 8.524997)
"8FVC9G8F+6X"

julia> # Encode a location using five digits of additional refinement:
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

© 2023 Klaus Crusius

[gha-img]: https://github.com/KlausC/OpenLocationCode.jl/actions/workflows/CI.yml/badge.svg?branch=main
[gha-url]: https://github.com/KlausC/OpenLocationCode.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/KlausC/OpenLocationCode.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/KlausC/OpenLocationCode.jl
