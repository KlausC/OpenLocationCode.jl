
export distance, area

const GEO_RADIUS = 6378.388e3 # meters

"""
    distance(lat1, lon1, lat2, lon2)

Distance between two points on a sphere with radius `r` in meters.
The point are given in degrees.
Uses haversine algorithm.
"""
function distance(lat1::Real, lon1::Real, lat2::Real, lon2::Real, r::Real=GEO_RADIUS)
    distance(promote(lat1, lon1, lat2, lon2, r)...)
end

function distance(lat1::T, lon1::T, lat2::T, lon2::T, r::T=6378.388e3) where T<:Real
    dlat = lat1 - lat2
    dlon = lon1 - lon2
    a = sind(dlat / 2) ^ 2 + sind(dlon / 2) ^ 2 * cosd(lat1) * cosd(lat2)
    dist =  asin(sqrt(a)) * r * 2
    dist
end

function distance(ca::CodeArea, cb::CodeArea)
    lata, lona = latitude_center(ca), longitude_center(ca)
    latb, lonb = latitude_center(cb), longitude_center(cb)
    distance(lata, lona, latb, lonb)
end

function distance(a::AbstractString, b::AbstractString)
    distance(decode(a), decode(b))
end

function area(ca::CodeArea)
    lat1 = latitude_high(ca)
    lat0 = latitude_low(ca)
    lod = longitude_high(ca) - longitude_low(ca)
    r = GEO_RADIUS
    area = (sind(lat1) - sind(lat0)) * lod / 180 * π * r^2
    area
end

area(ca::AbstractString) = area(decode(ca))
