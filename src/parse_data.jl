module Data

export network, tsunami_data, tsunami_extrema, tsunami_resolution, people, shelters;

import OpenStreetMapX;
import CSV;
import Makie;
import ArchGDAL;
const OSMX = OpenStreetMapX;
const AG = ArchGDAL;

include("utilities.jl");
using .Utils;

struct Person
    pos::Point
    evac::Bool
end

network(filename) = OSMX.get_map_data(filename, use_cache = false, trim_to_connected_graph = true);

# Get all of the tsunami inundation datasets as a dictionary;
# key is filename besides extension, value is data
datasets = Dict();
geotransform, offset, missing_data_val = nothing, nothing, nothing;

function tsunami_data(initial_time, half_mins, network, dir_name)
    ext = ".asc";
    for filename in filter(x -> endswith(x, ext), readdir(dir_name))
        datasets[filename[1:end-length(ext)]] = AG.readraster(dir_name * "/" * filename);
    end
    datasets[string(initial_time)] = datasets["30"];

    dataset = datasets[string(initial_time)];
    source = AG.importWKT(AG.getproj(dataset)); # For the tsunami asc files, the projection is UTM
    lla = AG.importEPSG(4326); # Lat-long
    global geotransform = AG.getgeotransform(dataset);

    # Determine how far off the tsunami dataset origin is from the road network origin, as ENU/UTM
    road_center = OSMX.center(network.bounds); # Lat-long of center of the road network, i.e., (0, 0) ENU
    network_point = AG.createpoint(road_center.lat, road_center.lon);
    AG.createcoordtrans(lla, source) do transform
        AG.transform!(network_point, transform)
    end
    # Assume everything besides the road network is in UTM
    global offset = (AG.getx(network_point, 0), AG.gety(network_point, 0));

    global missing_data_val = AG.getnodatavalue(AG.getband(datasets[string(initial_time)], 1));
    # Determine tsunami coordinates in ENU with the same origin as the road
    # network. UTM and ENU are both in meters, so no additional conversion needed
    x, y = begin
        """
        Get a list of possible values for a given dimension corresponding to x or y.
        """
        function tsunami_coord(dim, i)
            floats = Float64.(1:dim(datasets[string(initial_time)]));
            transform_order = i == 1 ? (floats, 1.0) : (1.0, floats);
            map(x -> x[i], AG.applygeotransform.(tuple(geotransform), transform_order...)) .- offset[i]
        end
        tsunami_coord(AG.width, 1), reverse(tsunami_coord(AG.height, 2))
    end
    # Replace the NaN value in the raw dataset (usually -9999) with NaN. Can't use
    # `missing` because the `heatmap` function doesn't know how to handle it. Also
    # flip the coordinates so it's not upside down (necessary because the asc files
    # store coordinates with an origin of top left, not bottom left)
    # Use half-mins because that's the interval of the tsunami data
    z = Makie.@lift(map(x -> x == missing_data_val ? NaN : x, reverse(datasets[string($half_mins * 30)][:, :, 1], dims=2)));
    x, y, z, offset
end

function tsunami_resolution()::NTuple{2,Float64}
    geotransform[2], -geotransform[6]
end

function tsunami_extrema()::NTuple{2,Float64}
    # Find minimum and maximum z values of the tsunami for heatmap normalization.
    # This can probably be done with `extrema()`, but it's nested so it's
    # complicated (and it's only computed once, unlike most calculations here)
    minᶻ = minimum(x -> minimum(filter(z -> z != missing_data_val, x[:, :, 1])), values(datasets));
    maxᶻ = maximum(x -> maximum(filter(z -> z != missing_data_val, x[:, :, 1])), values(datasets));
    minᶻ, maxᶻ
end

function people(filename)
    map(CSV.File(filename; normalizenames = true, types = [Int, Float64, Float64, Float64, Bool])) do p
        Person((p.X - offset[1], p.Y - offset[2]), p.Attribute_2)
    end
end

function shelters(network, filename)
    shelter_locs = CSV.File(filename; normalizenames = true, types = [Int, Float64, Float64]); # Shelter locations
    [begin
        node_id = OSMX.nearest_node(network, OSMX.ENU(shelter.x - offset[1], shelter.y - offset[2]));
        (node_id, enu_to_tuple(network.nodes[node_id]))
    end
    for shelter in shelter_locs]
end

function refresh_data(initial_time, half_mins, network_loc, tsunami_loc, people_loc, shelters_loc)
    new_network = network(network_loc);
    new_network,
    tsunami_data(initial_time, half_mins, new_network, tsunami_loc),
    tsunami_extrema(),
    people(people_loc),
    shelters(new_network, shelters_loc)
end

end # module
