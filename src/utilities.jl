module Utils

import OpenStreetMapX;
const OSMX = OpenStreetMapX;

export s, ft, mins, hr, mi, mph, m;
export Point, enu_to_tuple, dist;

# Unit conversions
# These take advantage of Julia's 'numeric literal coefficients' to appear like regular units
const s = 1.0; # Number of frames per second
const ft = 0.3048; # ft to m
const mins = 60s; # mins to frames
const hr = 60mins; # hrs to frames
const mi = 5280ft; # miles to m
const mph = 1mi / 1hr; # mph to m/frame
const m = 1.0; # Exactly 1, no conversion (included for completeness)

const Point = NTuple{2,Float64};

"""
Convert an ENU from the `OpenStreetMapX` library to a 2D tuple.
"""
enu_to_tuple(p::OSMX.ENU)::Point = (p.east, p.north);

"""
Distance between two points.
"""
dist(a::Point, b::Point)::Float64 = OSMX.distance(OSMX.ENU(a...), OSMX.ENU(b...));

"""
Determine the angle between two points.
"""
angle(a::Point, b::Point)::Float64 = atan(b[2] - a[2], b[1] - a[1]);

end # module
