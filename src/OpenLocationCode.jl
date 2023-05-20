module OpenLocationCode

export is_valid, is_short, is_full, encode, decode, recover_nearest, shorten
export CodeArea, latitude_low, longitude_low, latitude_high, longitude_high
export latitude_center, longitude_center, latitude_precision, longitude_precision

"A separator used to break the code into two parts to aid memorability."
const SEPARATOR_ = '+'

"The number of characters to place before the separator."
const SEPARATOR_POSITION_ = 8

# The character used to pad codes.
const PADDING_CHARACTER_ = '0'

# The character set used to encode the values.
const CODE_ALPHABET_ = "23456789CFGHJMPQRVWX"

# The base to use to convert numbers to/from.
const ENCODING_BASE_ = length(CODE_ALPHABET_) # 20

# The maximum value for latitude in degrees.
const LATITUDE_MAX_ = 90

# The maximum value for longitude in degrees.
const LONGITUDE_MAX_ = 180

# The max number of digits to process in a plus code.
const MAX_DIGIT_COUNT_ = 15

"""
Maximum code length using lat/lng pair encoding. The area of such a
code is approximately 13x13 meters (at the equator), and should be suitable
for identifying buildings. This excludes prefix and separator characters.
"""
const PAIR_CODE_LENGTH_ = 10

# First place value of the pairs (if the last pair value is 1).
const PAIR_FIRST_PLACE_VALUE_ = ENCODING_BASE_^(PAIR_CODE_LENGTH_ ÷ 2 - 1)

# Inverse of the precision of the pair section of the code.
const PAIR_PRECISION_ = ENCODING_BASE_^3

# The resolution values in degrees for each position in the lat/lng pair
# encoding. These give the place value of each position, and therefore the
# dimensions of the resulting area.
const PAIR_RESOLUTIONS_ = [20.0, 1.0, 0.05, 0.0025, 0.000125]

# Number of digits in the grid precision part of the code.
const GRID_CODE_LENGTH_ = MAX_DIGIT_COUNT_ - PAIR_CODE_LENGTH_

# Number of columns in the grid refinement method.
const GRID_COLUMNS_ = 4

# Number of rows in the grid refinement method.
const GRID_ROWS_ = 5

# First place value of the latitude grid (if the last place is 1).
const GRID_LAT_FIRST_PLACE_VALUE_ = GRID_ROWS_^(GRID_CODE_LENGTH_ - 1)

# First place value of the longitude grid (if the last place is 1).
const GRID_LNG_FIRST_PLACE_VALUE_ = GRID_COLUMNS_^(GRID_CODE_LENGTH_ - 1)

# inverse of best precision of grid part
const GRID_LAT_PRECISION = GRID_LAT_FIRST_PLACE_VALUE_ * GRID_ROWS_
const GRID_LNG_PRECISION = GRID_LNG_FIRST_PLACE_VALUE_ * GRID_COLUMNS_

# Multiply latitude by this much to make it a multiple of the finest
# precision.
const FINAL_LAT_PRECISION_ = PAIR_PRECISION_ * GRID_LAT_PRECISION

# Multiply longitude by this much to make it a multiple of the finest
# precision.
const FINAL_LNG_PRECISION_ = PAIR_PRECISION_ * GRID_LNG_PRECISION

# Minimum length of a code that can be shortened.
const MIN_TRIMMABLE_CODE_LEN_ = 6

const GRID_SIZE_DEGREES_ = 0.000125

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

min_lat(ca, delta) = min(latitude_low(ca) + delta, LATITUDE_MAX_ )
min_lon(ca, delta) = min(longitude_low(ca) + delta, LONGITUDE_MAX_)

function Base.isapprox(ca::CodeArea, cb::CodeArea)
    ca.codelength == cb.codelength && ca.latlo ≈ cb.latlo && ca.longlo ≈ cb.longlo
end

"""
    is_valid(code::AbstractString)

Determine if a code is valid.
To be valid, all characters must be from the Open Location Code character
set with at most one separator. The separator can be in any even-numbered
position up to the eighth digit.
"""
function is_valid(code::AbstractString)

    # The separator is required.
    seps = first.(findall(SEPARATOR_ * "", code))
    length(seps) != 1 && return false
    sep = seps[begin]
    sep - 1 > SEPARATOR_POSITION_ && return false
    # Is it the only character?
    length(code) <= 1 && return false
    # Is it in an illegal position?
    sep % 2 == 0 && return false
    # We can have an even number of padding characters before the separator,
    # but then it must be the final character.
    pads = first.(findall(PADDING_CHARACTER_ * "", code))
    if !isempty(pads)
        pad = pads[begin]
        # Short codes cannot have padding
        sep < SEPARATOR_POSITION_ && return false
        # Not allowed to start with them!
        pad == 1 && return false

        # There can only be one group and it must have even length.
        (length(pads) % 2 != 0 || pads[end] - pad + 1 != length(pads)) && return false
        # If the code is long enough to end with a separator, make sure it does.
        return code[end] == SEPARATOR_
    end
    # If there are characters after the separator, make sure there isn't just
    # one of them (not legal).
    length(code) - sep == 1 && return false
    # Check the code contains only valid characters.
    sepPad = SEPARATOR_ * PADDING_CHARACTER_
    return all(ch -> uppercase(ch) in CODE_ALPHABET_ || ch in sepPad, code)
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
    sep = findfirst(SEPARATOR_, code)
    return sep <= SEPARATOR_POSITION_
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
    firstLatValue = (findfirst(uppercase(code[1]), CODE_ALPHABET_) - 1) * ENCODING_BASE_
    if firstLatValue >= LATITUDE_MAX_ * 2
        # The code would decode to a latitude of >= 90 degrees.
        return false
    end
    # Work out what the first longitude character indicates for longitude.
    firstLngValue = (findfirst(uppercase(code[2]), CODE_ALPHABET_) - 1) * ENCODING_BASE_
    if firstLngValue >= LONGITUDE_MAX_ * 2
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
10 characters, returning a code of approximately 13.5x13.5 meters. Longer
codes represent smaller areas, but lengths > 14 are sub-centimetre and so
11 or 12 are probably the limit of useful codes.
Args:
  latitude: A latitude in signed decimal degrees. Will be clipped to the
      range -90 to 90.
  longitude: A longitude in signed decimal degrees. Will be normalised to
      the range -180 to 180.
  codelength: The number of significant digits in the output code, not
      including any separator characters.
"""
function encode(latitude::Real, longitude::Real, codelength=PAIR_CODE_LENGTH_)
    encode(float.(promote(latitude, longitude))..., codelength)
end

function encode(latitude::T, longitude::T, codelength::Int=PAIR_CODE_LENGTH_) where {T<:AbstractFloat}
    if codelength < 2 || (codelength < PAIR_CODE_LENGTH_ && codelength % 2 == 1)
        throw(ArgumentError("Invalid Open Location Code length - $codelength"))
    end
    codelength = min(codelength, MAX_DIGIT_COUNT_)
    # Ensure that latitude and longitude are valid.
    latitude = clipLatitude(latitude)
    longitude = normalizeLongitude(longitude)
    # Latitude 90 needs to be adjusted to be just less, so the returned code
    # can also be decoded.
    if latitude == LATITUDE_MAX_
        latitude = latitude - latitude_precision(codelength)
    end
    code = ""

    # Compute the code.
    # This approach converts each value to an integer after multiplying it by
    # the final precision. This allows us to use only integer operations, so
    # avoiding any accumulation of floating point representation errors.

    # Multiply values by their precision and convert to positive.
    # Force to integers so the division operations will have integer results.
    latVal = valpairgrid(latitude, FINAL_LAT_PRECISION_, LATITUDE_MAX_)
    lngVal = valpairgrid(longitude, FINAL_LNG_PRECISION_, LONGITUDE_MAX_)

    # Compute the grid part of the code if necessary.
    if codelength > PAIR_CODE_LENGTH_
        for _ in 1:GRID_CODE_LENGTH_
            latDigit = latVal % GRID_ROWS_
            lngDigit = lngVal % GRID_COLUMNS_
            ndx = latDigit * GRID_COLUMNS_ + lngDigit
            code = CODE_ALPHABET_[ndx+1] * code
            latVal ÷= GRID_ROWS_
            lngVal ÷= GRID_COLUMNS_
        end
    else
        latVal ÷= ^(GRID_ROWS_, GRID_CODE_LENGTH_)
        lngVal ÷= ^(GRID_COLUMNS_, GRID_CODE_LENGTH_)
    end
    # Compute the pair section of the code.
    for _ in 1:PAIR_CODE_LENGTH_÷2
        code = CODE_ALPHABET_[lngVal%ENCODING_BASE_+1] * code
        code = CODE_ALPHABET_[latVal%ENCODING_BASE_+1] * code
        latVal ÷= ENCODING_BASE_
        lngVal ÷= ENCODING_BASE_
    end

    # Add the separator character.
    code = code[1:SEPARATOR_POSITION_] * SEPARATOR_ * code[SEPARATOR_POSITION_+1:end]

    # If we don't need to pad the code, return the requested section.
    if codelength >= SEPARATOR_POSITION_
        return code[1:codelength+1]
    end

    # Pad and return the code.
    return code[1:codelength] * '0'^(SEPARATOR_POSITION_ - codelength) * SEPARATOR_
end

"""
    valpairgrid(lat, gridprecision)

Return integer value in multiples of gridprecision
"""
function valpairgrid(latlo::AbstractFloat, prec::Integer, max::Integer)
    unsafe_trunc(Int64, floor(latlo * prec + 2.0^-23)) + Int64(prec) * max
end

"""
    decode(code)

Decode an Open Location Code into the location coordinates.
Returns a CodeArea object that includes the coordinates of the bounding
box - the lower left, center and upper right.
Args:
  code: The Open Location Code to decode.
Returns:
  A CodeArea object that provides the latitude and longitude of two of the
  corners of the area, the center, and the length of the original code.
"""
function decode(code::AbstractString)

    if !is_full(code)
        throw(ArgumentError("Passed Open Location Code is not a valid full code - $code"))
    end
    # Strip out separator character (we've already established the code is
    # valid so the maximum is one), and padding characters. Convert to upper
    # case and constrain to the maximum number of digits.
    m = match(r"([^0+]*)[0]*[+]?(.*)", code)
    code = uppercase(m[1] * m[2])
    if length(code) > MAX_DIGIT_COUNT_
        code = code[1:MAX_DIGIT_COUNT_]
    end
    # Initialise the values for each section. We work them out as integers and
    # convert them to floats at the end.
    normalLat = -LATITUDE_MAX_ * PAIR_PRECISION_
    normalLng = -LONGITUDE_MAX_ * PAIR_PRECISION_
    gridLat = 0
    gridLng = 0
    # How many digits do we have to process?
    digits = min(length(code), PAIR_CODE_LENGTH_)
    # Define the place value for the most significant pair.
    pv = PAIR_FIRST_PLACE_VALUE_
    # Decode the paired digits.
    for i in 1:2:digits
        normalLat += (findfirst(code[i], CODE_ALPHABET_) - 1) * pv
        normalLng += (findfirst(code[i+1], CODE_ALPHABET_) - 1) * pv
        if i < digits - 2
            pv ÷= ENCODING_BASE_
        end
    end
    # Convert the place value to a float in degrees.
    latPrecision = pv / PAIR_PRECISION_
    lngPrecision = pv / PAIR_PRECISION_
    # Process any extra precision digits.
    if length(code) > PAIR_CODE_LENGTH_
        # Initialise the place values for the grid.
        rowpv = GRID_LAT_FIRST_PLACE_VALUE_
        colpv = GRID_LNG_FIRST_PLACE_VALUE_
        # How many digits do we have to process?
        digits = min(length(code), MAX_DIGIT_COUNT_)
        for i in PAIR_CODE_LENGTH_:digits-1
            digitVal = findfirst(code[i+1], CODE_ALPHABET_) - 1
            row = digitVal ÷ GRID_COLUMNS_
            col = digitVal % GRID_COLUMNS_
            gridLat += row * rowpv
            gridLng += col * colpv
            if i < digits - 1
                rowpv ÷= GRID_ROWS_
                colpv ÷= GRID_COLUMNS_
            end
        end
        # Adjust the precisions from the integer values to degrees.
        latPrecision = rowpv / FINAL_LAT_PRECISION_
        lngPrecision = colpv / FINAL_LNG_PRECISION_
    end
    # Merge the values from the normal and extra precision parts of the code.
    lat = normalLat / PAIR_PRECISION_ + gridLat / FINAL_LAT_PRECISION_
    lng = normalLng / PAIR_PRECISION_ + gridLng / FINAL_LNG_PRECISION_
    # Multiply values by 1e14, round and then divide. This reduces errors due
    # to floating point precision.
    digits = 14
    base = 10
    CodeArea(round(lat; digits, base),
        round(lng; digits, base),
        min(length(code), MAX_DIGIT_COUNT_))
end

"""
    recover_nearest(code, latitude, longitude)

Recover the nearest matching code to a specified location.
Given a short code of between four and seven characters, this recovers
the nearest matching full code to the specified location.
Args:
  code: A valid OLC character sequence.
  referenceLatitude: The latitude (in signed decimal degrees) to use to
      find the nearest matching full code.
  referenceLongitude: The longitude (in signed decimal degrees) to use
      to find the nearest matching full code.
Returns:
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
    paddingLength = SEPARATOR_POSITION_ - findfirst(SEPARATOR_, code) + 1
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
       latcenter - resolution >= -LATITUDE_MAX_
        # If the proposed code is more than half a cell north of the reference location,
        # it's too far, and the best match will be one cell south.
        latcenter -= resolution
    elseif referenceLatitude - halfResolution > latcenter &&
           latcenter + resolution <= LATITUDE_MAX_
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
Args:
  code: A full, valid code to shorten.
  latitude: A latitude, in signed decimal degrees, to use as the reference
      point.
  longitude: A longitude, in signed decimal degrees, to use as the reference
      point.
Returns:
  Either the original code, if the reference location was not close enough,
  or the .
"""
function shorten(code, latitude, longitude)

    if !is_full(code)
        throw(ArgumentError("Passed code is not valid and full: $code"))
    end
    if findfirst(PADDING_CHARACTER_, code) !== nothing
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
    for i in length(PAIR_RESOLUTIONS_)-2:-1:0
        # Check if we're close enough to shorten. The range must be less than 1/2
        # the resolution to shorten at all, and we want to allow some safety, so
        # use 0.3 instead of 0.5 as a multiplier.
        if coderange < (PAIR_RESOLUTIONS_[i+1] * 0.3)
            # Trim it.
            return code[(i+1)*2+1:end]
        end
    end
    return code
end

"""
    clipLatitude(latitude)

Clip a latitude into the range -90 to 90.
Args:
  latitude: A latitude in signed decimal degrees.
"""
function clipLatitude(latitude::Real)
    return min(LATITUDE_MAX_, max(-LATITUDE_MAX_, latitude))
end

"""
    computeLatitudePrecision(codelength)

Compute the latitude precision value for a given code length. Lengths <=
10 have the same precision for latitude and longitude, but lengths > 10
have different precisions due to the grid method having fewer columns than
rows.
"""
latitude_precision(codelength) = compute_precision(codelength, GRID_ROWS_)
longitude_precision(codelength) = compute_precision(codelength, GRID_COLUMNS_)

function compute_precision(codelength, grid)
    if codelength <= PAIR_CODE_LENGTH_
        ENCODING_BASE_ / (ENCODING_BASE_^(codelength ÷ 2 - 1))
    else
        1 / (ENCODING_BASE_^3 * grid^(codelength - PAIR_CODE_LENGTH_))
    end
end

"""
    normalizeLongitude(longitude)

Normalize a longitude into the range -180 to 180, not including 180.
Args:
  longitude: A longitude in signed decimal degrees.
"""
function normalizeLongitude(longitude::Real)
    while longitude < -LONGITUDE_MAX_
        longitude = longitude + LONGITUDE_MAX_ * 2
    end
    while longitude >= LONGITUDE_MAX_
        longitude = longitude - LONGITUDE_MAX_ * 2
    end
    return longitude
end

end
