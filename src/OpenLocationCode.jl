module OpenLocationCode

export is_valid, is_short, is_full, encode, decode, recover_nearest, shorten
export CodeArea, latitude_low, longitude_low, latitude_high, longitude_high
export latitude_center, longitude_center, latitude_precision, longitude_precision

"A separator used to break the code into two parts to aid memorability."
const SEPARATOR = '+'

"The number of characters to place before the separator."
const SEPARATOR_POSITION = 8

# The character used to pad codes.
const PADDING_CHARACTER = '0'

# The character set used to encode the values.
const CODE_DIGITS = [UInt8(c) for c in "23456789CFGHJMPQRVWX"]

# The base to use to convert numbers to/from.
const ENCODING_BASE = length(CODE_DIGITS) # 20

# The maximum value for latitude in degrees.
const LATITUDE_MAX = 90

# The maximum value for longitude in degrees.
const LONGITUDE_MAX = 180

# The max number of digits to process in a plus code.
const MAX_DIGIT_COUNT = 15

"""
Maximum code length using lat/lng pair encoding. The area of such a
code is approximately 13x13 meters (at the equator), and should be suitable
for identifying buildings. This excludes prefix and separator characters.
"""
const PAIR_CODE_LENGTH = 10

# First place value of the pairs (if the last pair value is 1).
const PAIR_FIRST_PLACE_VALUE = ENCODING_BASE^(PAIR_CODE_LENGTH ÷ 2 - 1)

# Inverse of the precision of the pair section of the code.
const PAIR_PRECISION = ENCODING_BASE^3

# The resolution values in degrees for each position in the lat/lng pair
# encoding. These give the place value of each position, and therefore the
# dimensions of the resulting area.
const PAIR_RESOLUTIONS = [20.0, 1.0, 0.05, 0.0025, 0.000125]

# Number of digits in the grid precision part of the code.
const GRID_CODE_LENGTH = MAX_DIGIT_COUNT - PAIR_CODE_LENGTH

# Number of columns in the grid refinement method.
const GRID_COLUMNS = 4

# Number of rows in the grid refinement method.
const GRID_ROWS = 5

# First place value of the latitude grid (if the last place is 1).
const GRID_LAT_FIRST_PLACE_VALUE = GRID_ROWS^(GRID_CODE_LENGTH - 1)

# First place value of the longitude grid (if the last place is 1).
const GRID_LNG_FIRST_PLACE_VALUE = GRID_COLUMNS^(GRID_CODE_LENGTH - 1)

# inverse of best precision of grid part
const GRID_LATPRECISION = GRID_LAT_FIRST_PLACE_VALUE * GRID_ROWS
const GRID_LNGPRECISION = GRID_LNG_FIRST_PLACE_VALUE * GRID_COLUMNS

# Multiply latitude by this much to make it a multiple of the finest
# precision.
const FINAL_LAT_PRECISION = PAIR_PRECISION * GRID_LATPRECISION

# Multiply longitude by this much to make it a multiple of the finest
# precision.
const FINAL_LNG_PRECISION = PAIR_PRECISION * GRID_LNGPRECISION

# Minimum length of a code that can be shortened.
const MIN_TRIMMABLE_CODE_LEN = 6

const GRID_SIZE_DEGREES = 0.000125

"""
    CodeArea(latLo, longLo, codelength)

Coordinates of a decoded Open Location Code.
The coordinates include the latitude and longitude of the lower left
and the code length, which defines the extension of the area.

Attributes:
  latitude_low: The latitude of the SW corner in degrees.
  longitude_low: The longitude of the SW corner in degrees.
  latitude_high: The latitude of the NE corner in degrees.
  longitude_high: The longitude of the NE corner in degrees.
  latitude_center: The latitude of the center in degrees.
  longitude_center: The longitude of the center in degrees.
  code_length: The number of significant characters that were in the code.
"""
struct CodeArea{T<:Real}
    latlo::T
    longlo::T
    codelength::Int

    function CodeArea(latlo, longlo, codelength)
        latlo, longlo = promote(latlo, longlo)
        new{typeof(latlo)}(latlo, longlo, codelength)
    end
end

latitude_low(ca::CodeArea) = ca.latlo
longitude_low(ca::CodeArea) = ca.longlo
latitude_high(ca::CodeArea) = min_lat(ca, latitude_precision(ca))
longitude_high(ca::CodeArea) = min_lon(ca, longitude_precision(ca))
latitude_center(ca::CodeArea) = min_lat(ca, latitude_precision(ca) / 2)
longitude_center(ca::CodeArea) = min_lon(ca, longitude_precision(ca) / 2)
latitude_precision(ca::CodeArea) = latitude_precision(ca.codelength)
longitude_precision(ca::CodeArea) = longitude_precision(ca.codelength)

min_lat(ca, delta) = min(latitude_low(ca) + delta, LATITUDE_MAX )
min_lon(ca, delta) = min(longitude_low(ca) + delta, LONGITUDE_MAX)

function Base.isapprox(ca::CodeArea, cb::CodeArea)
    ca.codelength == cb.codelength && ca.latlo ≈ cb.latlo && ca.longlo ≈ cb.longlo
end

function Base.show(io::IO, ca::CodeArea)
    print(io, "CodeArea(")
    print(io, latitude_low(ca), "+", latitude_precision(ca), ", ")
    print(io, longitude_low(ca), "+", longitude_precision(ca), ", ")
    print(io, ca.codelength, ")")
end

"""
    is_valid(code::AbstractString)

Determine if a code is valid.
To be valid, all characters must be from the Open Location Code character
set with at most one separator. The separator can be in any even-numbered
position up to the eighth digit.
"""
function is_valid(code::AbstractString)

    # The separator is required exactly once.
    sep = findfirst(SEPARATOR, code)
    sep !== nothing && findlast(SEPARATOR, code) == sep || return false
    # Not more than 8 characters before separator.
    sep > SEPARATOR_POSITION + 1 && return false
    # Is it the only character?
    length(code) <= 1 && return false
    # Is it in an illegal position?
    sep % 2 == 0 && return false
    # We can have an even number of padding characters before the separator,
    # but then it must be the final character.
    pad = findfirst(PADDING_CHARACTER, code)
    if pad !== nothing
        # Short codes cannot have padding
        sep < SEPARATOR_POSITION && return false
        # Not allowed to start with them or start at even poaition
        pad == 1 || iseven(pad) && return false

        # There can only be one group and it must have even length.
        pad2 = findlast(PADDING_CHARACTER, code)
        (pad2 - pad) % 2 == 0 && return false
        all(isequal(PADDING_CHARACTER), view(code, pad:pad2)) || return false
        # If the code is long enough to end with a separator, make sure it does.
        return code[end] == SEPARATOR
    end
    # If there are characters after the separator, make sure there isn't just
    # one of them (not legal).
    length(code) - sep == 1 && return false
    # Check the code contains only valid characters.
    return all(ch -> UInt8(uppercase(ch)) in CODE_DIGITS || ch == SEPARATOR || ch == PADDING_CHARACTER, code)
end

"""
    is_short(code)

Determine if a code is a valid short code.
A short Open Location Code is a sequence created by removing four or more
digits from an Open Location Code. It must include a separator
character.
"""
function is_short(code::AbstractString)
    # Check it's valid.
    is_valid(code) || return false
    # If there are less characters than expected before the SEPARATOR.
    sep = findfirst(SEPARATOR, code)
    return sep <= SEPARATOR_POSITION
end

"""
    is_full(code)

Determine if a code is a valid full Open Location Code.
Not all possible combinations of Open Location Code characters decode to
valid latitude and longitude values. This checks that a code is valid
and also that the latitude and longitude values are legal. If the prefix
character is present, it must be the first character. If the separator
character is present, it must be after four characters.
"""
function is_full(code::AbstractString)
    is_valid(code) || return false
    # If it's short, it's not full
    is_short(code) && return false
    # Work out what the first latitude character indicates for latitude.
    clat = UInt8(uppercase(code[1]))
    firstLatValue = (findfirst(==(clat), CODE_DIGITS) - 1) * ENCODING_BASE
    if firstLatValue >= LATITUDE_MAX * 2
        # The code would decode to a latitude of >= 90 degrees.
        return false
    end
    # Work out what the first longitude character indicates for longitude.
    clng = UInt8(uppercase(code[2]))
    firstLngValue = (findfirst(==(clng), CODE_DIGITS) - 1) * ENCODING_BASE
    if firstLngValue >= LONGITUDE_MAX * 2
        # The code would decode to a longitude of >= 180 degrees.
        return false
    end
    return true
end

"""
    encode(latitude, longitude[, codelength])

Encode a location into an Open Location Code.
Produces a code of the specified length, or the default length if no length
is provided.

The length determines the accuracy of the code. The default length is
10 characters, returning a code of approximately 13.9x13.9 meters. Longer
codes represent smaller areas, but lengths > 15 are sub-centimetre and so
11 or 12 are probably the limit of useful codes.

# Arguments
- `latitude`: A latitude in signed degrees. Will be clipped to the range -90 to 90.
- `longitude`: A longitude in signed degrees. Will be normalised to
    the range -180 to 180.
- `codelength`: The number of significant digits in the output code, not
    including any separator characters.

# Examples:
```julia
julia> encode(50.173168, 8.338086, 11)
"9F2C58FQ+768"
```
"""
function encode(latitude::Real, longitude::Real, codelength=PAIR_CODE_LENGTH)
    encode(float.(promote(latitude, longitude))..., codelength)
end

function encode(latitude::T, longitude::T, codelength::Int=PAIR_CODE_LENGTH) where {T<:AbstractFloat}
    if codelength < 2 || (codelength < PAIR_CODE_LENGTH && codelength % 2 == 1)
        throw(ArgumentError("Invalid Open Location Code length - $codelength"))
    end
    codelength = min(codelength, MAX_DIGIT_COUNT)
    # Ensure that latitude and longitude are valid.
    latitude = clipLatitude(latitude)
    longitude = normalizeLongitude(longitude)
    # Latitude 90 needs to be adjusted to be just less, so the returned code
    # can also be decoded.
    if latitude == LATITUDE_MAX
        latitude = latitude - latitude_precision(codelength)
    end
    bcode = UInt8[]

    # Compute the code.
    # This approach converts each value to an integer after multiplying it by
    # the final precision. This allows us to use only integer operations, so
    # avoiding any accumulation of floating point representation errors.

    # Multiply values by their precision and convert to positive.
    # Force to integers so the division operations will have integer results.
    latVal = valpairgrid(latitude, FINAL_LAT_PRECISION, LATITUDE_MAX)
    lngVal = valpairgrid(longitude, FINAL_LNG_PRECISION, LONGITUDE_MAX)

    # Compute the grid part of the code if necessary.
    if codelength > PAIR_CODE_LENGTH
        for _ in 1:GRID_CODE_LENGTH
            latDigit = latVal % GRID_ROWS
            lngDigit = lngVal % GRID_COLUMNS
            ndx = latDigit * GRID_COLUMNS + lngDigit
            push!(bcode, CODE_DIGITS[ndx+1])
            latVal ÷= GRID_ROWS
            lngVal ÷= GRID_COLUMNS
        end
    else
        latVal ÷= ^(GRID_ROWS, GRID_CODE_LENGTH)
        lngVal ÷= ^(GRID_COLUMNS, GRID_CODE_LENGTH)
    end
    # Compute the pair section of the code.
    for _ in 1:PAIR_CODE_LENGTH÷2
        push!(bcode, CODE_DIGITS[lngVal%ENCODING_BASE+1])
        push!(bcode, CODE_DIGITS[latVal%ENCODING_BASE+1])
        latVal ÷= ENCODING_BASE
        lngVal ÷= ENCODING_BASE
    end

    reverse!(bcode)
    # If we don't need to pad the code, return the requested section.
    if codelength >= SEPARATOR_POSITION
        resize!(bcode, codelength)
        insert!(bcode, SEPARATOR_POSITION + 1, SEPARATOR)
        return String(bcode)
    end

    # Pad and return the code.
    resize!(bcode, codelength)
    for i = codelength+1:SEPARATOR_POSITION
        push!(bcode, PADDING_CHARACTER)
    end
    push!(bcode, SEPARATOR)
    return String(bcode)
end

"""
    valpairgrid(lat, precision, max)

Return integer value in multiples of `1/precision`
"""
function valpairgrid(latlo::AbstractFloat, prec::Integer, max::Integer)
    unsafe_trunc(Int64, floor(nextfloat(latlo) * prec)) + Int64(prec) * max
end

"""
    decode(code)

Decode an Open Location Code into the location coordinates.

# Arguments:
- `code`: The Open Location Code to decode.

# Returns:
  A `CodeArea` object that provides the latitude and longitude of two of the
  corners of the area, the center, and the length of the original code.
"""
function decode(code::AbstractString)

    if !is_full(code)
        throw(ArgumentError("Passed Open Location Code is not a valid full code - $code"))
    end
    # Strip out separator character (we've already established the code is
    # valid so the maximum is one), and padding characters. Convert to upper
    # case and constrain to the maximum number of digits.
    bcode = UInt8[]
    for c in code
        c == PADDING_CHARACTER && break
        c == SEPARATOR && continue
        push!(bcode, uppercase(c))
    end
    if length(bcode) > MAX_DIGIT_COUNT
        resize!(bcode, MAX_DIGIT_COUNT)
    end
    codelength = length(bcode)
    # Initialise the values for each section. We work them out as integers and
    # convert them to floats at the end.
    normalLat = -LATITUDE_MAX * PAIR_PRECISION
    normalLng = -LONGITUDE_MAX * PAIR_PRECISION
    gridLat = 0
    gridLng = 0
    # How many digits do we have to process?
    digits = min(codelength, PAIR_CODE_LENGTH)
    # Define the place value for the most significant pair.
    pv = PAIR_FIRST_PLACE_VALUE
    # Decode the paired digits.
    for i in 1:2:digits
        normalLat += (findfirst(==(bcode[i]), CODE_DIGITS) - 1) * pv
        normalLng += (findfirst(==(bcode[i+1]), CODE_DIGITS) - 1) * pv
        if i < digits - 2
            pv ÷= ENCODING_BASE
        end
    end
    # Convert the place value to a float in degrees.
    latPrecision = pv / PAIR_PRECISION
    lngPrecision = pv / PAIR_PRECISION
    # Process any extra precision digits.
    if codelength > PAIR_CODE_LENGTH
        # Initialise the place values for the grid.
        rowpv = GRID_LAT_FIRST_PLACE_VALUE
        colpv = GRID_LNG_FIRST_PLACE_VALUE
        # How many digits do we have to process?
        digits = min(codelength, MAX_DIGIT_COUNT)
        for i in PAIR_CODE_LENGTH:digits-1
            digitVal = findfirst(==(bcode[i+1]), CODE_DIGITS) - 1
            row = digitVal ÷ GRID_COLUMNS
            col = digitVal % GRID_COLUMNS
            gridLat += row * rowpv
            gridLng += col * colpv
            if i < digits - 1
                rowpv ÷= GRID_ROWS
                colpv ÷= GRID_COLUMNS
            end
        end
        # Adjust the precisions from the integer values to degrees.
        latPrecision = rowpv / FINAL_LAT_PRECISION
        lngPrecision = colpv / FINAL_LNG_PRECISION
    end
    # Merge the values from the normal and extra precision parts of the code.
    lat = normalLat / PAIR_PRECISION + gridLat / FINAL_LAT_PRECISION
    lng = normalLng / PAIR_PRECISION + gridLng / FINAL_LNG_PRECISION
    CodeArea(lat, lng, codelength)
end

"""
    recover_nearest(code, latitude, longitude)

Recover the nearest matching code to a specified location.
Given a short code of between four and seven characters, this recovers
the nearest matching full code to the specified location.

# Arguments:
- `code`: A valid OLC character sequence.
- `latitude`: The latitude (in signed degrees) to use to
      find the nearest matching full code.
- `longitude``: The longitude (in signed degrees) to use
      to find the nearest matching full code.

# Returns:
  The nearest full Open Location Code to the reference location that matches
  the short code. If the passed code was not a valid short code, but was a
  valid full code, it is returned with proper capitalization but otherwise
  unchanged.
"""
function recover_nearest(code::AbstractString, latitude::Real, longitude::Real)

    # if code is a valid full code, return it properly capitalized
    is_full(code) && return uppercase(code)
    if !is_short(code)
        throw(ArgumentError("Passed short code is not valid - $code"))
    end
    # Ensure that latitude and longitude are valid.
    referenceLatitude = clipLatitude(latitude)
    referenceLongitude = normalizeLongitude(longitude)
    # Clean up the passed code.
    code = uppercase(code)
    # Compute the number of digits we need to recover.
    paddingLength = SEPARATOR_POSITION - findfirst(SEPARATOR, code) + 1
    # The resolution (height and width) of the padded area in degrees.
    resolution = compute_precision(paddingLength, 1)
    # Distance from the center to an edge (in degrees).
    halfResolution = resolution / 2
    # Use the reference location to pad the supplied short code and decode it.
    codeArea = decode(encode(referenceLatitude, referenceLongitude)[1:paddingLength] * code)
    # How many degrees latitude is the code from the reference? If it is more
    # than half the resolution, we need to move it north or south but keep it
    # within -90 to 90 degrees.
    latcenter = latitude_center(codeArea)
    if referenceLatitude + halfResolution < latcenter &&
       latcenter - resolution >= -LATITUDE_MAX
        # If the proposed code is more than half a cell north of the reference location,
        # it's too far, and the best match will be one cell south.
        latcenter -= resolution
    elseif referenceLatitude - halfResolution > latcenter &&
           latcenter + resolution <= LATITUDE_MAX
        # If the proposed code is more than half a cell south of the reference location,
        # it's too far, and the best match will be one cell north.
        latcenter += resolution
    end
    # Adjust longitude if necessary.
    longcenter = longitude_center(codeArea)
    if referenceLongitude + halfResolution < longcenter
        longcenter -= resolution
    elseif referenceLongitude - halfResolution > longcenter
        longcenter += resolution
    end
    return encode(latcenter, longcenter, codeArea.codelength)
end

"""
    shorten(code, latitude, longitude)

Remove characters from the start of an OLC code.
This uses a reference location to determine how many initial characters
can be removed from the OLC code. The number of characters that can be
removed depends on the distance between the code center and the reference
location.

The minimum number of characters that will be removed is four. If more than
four characters can be removed, the additional characters will be replaced
with the padding character. At most eight characters will be removed.
The reference location must be within 50% of the maximum range. This ensures
that the shortened code will be able to be recovered using slightly different
locations.

# Arguments

- `code`: A full, valid code to shorten.
- `latitude`: A latitude, in signed degrees, to use as the reference point.
- `longitude`: A longitude, in signed degrees, to use as the reference point.

# Returns:
  Either the original code, if the reference location was not close enough,
  or the shortest code which can be used to recover from the reference location.
"""
function shorten(code, latitude, longitude)

    if !is_full(code)
        throw(ArgumentError("Passed code is not valid and full: $code"))
    end
    if findfirst(PADDING_CHARACTER, code) !== nothing
        throw(ArgumentError("Cannot shorten padded codes: $code"))
    end
    code = uppercase(code)
    codeArea = decode(code)
    # Ensure that latitude and longitude are valid.
    latitude = clipLatitude(latitude)
    longitude = normalizeLongitude(longitude)
    # How close are the latitude and longitude to the code center.
    coderange = max(abs(latitude_center(codeArea) - latitude),
        abs(longitude_center(codeArea) - longitude))
    for i in length(PAIR_RESOLUTIONS)-2:-1:0
        # Check if we're close enough to shorten. The range must be less than 1/2
        # the resolution to shorten at all, and we want to allow some safety, so
        # use 0.3 instead of 0.5 as a multiplier.
        if coderange < (PAIR_RESOLUTIONS[i+1] * 0.3)
            # Trim it.
            return code[(i+1)*2+1:end]
        end
    end
    return code
end

"""
    computeLatitudePrecision(codelength)

Compute the latitude precision value for a given code length. Lengths <=
10 have the same precision for latitude and longitude, but lengths > 10
have different precisions due to the grid method having fewer columns than
rows.
"""
latitude_precision(codelength) = compute_precision(codelength, GRID_ROWS)
longitude_precision(codelength) = compute_precision(codelength, GRID_COLUMNS)

function compute_precision(codelength, grid)
    if codelength <= PAIR_CODE_LENGTH
        ENCODING_BASE / (ENCODING_BASE^(codelength ÷ 2 - 1))
    else
        1 / (ENCODING_BASE^3 * grid^(codelength - PAIR_CODE_LENGTH))
    end
end

"""
    clipLatitude(latitude)

Clip a latitude into the range -90 to 90.
Args:
  latitude: A latitude in signed degrees.
"""
function clipLatitude(latitude::Real)
    return min(LATITUDE_MAX, max(-LATITUDE_MAX, latitude))
end

"""
    normalizeLongitude(longitude)

Normalize a longitude into the range -180 to 180, not including 180.
Args:
  longitude: A longitude in signed degrees.
"""
function normalizeLongitude(longitude::Real)
    while longitude < -LONGITUDE_MAX
        longitude = longitude + LONGITUDE_MAX * 2
    end
    while longitude >= LONGITUDE_MAX
        longitude = longitude - LONGITUDE_MAX * 2
    end
    return longitude
end

end
