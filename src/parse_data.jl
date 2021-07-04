module Data

export network, set_elevations!, get_slopes, get_landslides, tsunami_data, tsunami_extrema, tsunami_resolution, people, shelters;

import OpenStreetMapX;
import CSV;
import GLMakie;
import ArchGDAL;
const OSMX = OpenStreetMapX;
const AG = ArchGDAL;

include("utilities.jl");
using .Utils;

struct Person
    id::Int
    pos::Point
    evac::Bool
end

mutable struct Shelter
    id::Int
    node_id::Int
    pos::Point
    inundated::Bool
end

function network(filename)
    net = OSMX.get_map_data(filename, use_cache = false, trim_to_connected_graph = true);
    net, Dict((seg.node0, seg.node1) => seg.nodes for seg in OSMX.find_segments(net.nodes, net.roadways, net.intersections))
end

offset = nothing;

function read_asc(network, filepath)
    dataset = AG.readraster(filepath);

    source = AG.importWKT(AG.getproj(dataset)); # For the asc files, the projection is UTM
    lla = AG.importEPSG(4326); # Lat-long
    geotrans = AG.getgeotransform(dataset);

    # Determine how far off the dataset origin is from the road network origin, as ENU/UTM
    road_center = OSMX.center(network.bounds); # Lat-long of center of the road network, i.e., (0, 0) ENU
    network_point = AG.createpoint(road_center.lat, road_center.lon);
    AG.createcoordtrans(lla, source) do transform
        AG.transform!(network_point, transform)
    end
    # Assume people and shelters are the same coordinate system as elevation and landslide data
    global offset = (AG.getx(network_point, 0), AG.gety(network_point, 0));

    missing_val = AG.getnodatavalue(AG.getband(dataset, 1));
    # Determine coordinates in ENU with the same origin as the road network. UTM
    # and ENU are both in meters, so no additional conversion needed
    x_step, y_step = geotrans[2], -geotrans[6];
    origin = (geotrans[1] - offset[1] + x_step,
              geotrans[4] - offset[2] - y_step * size(dataset, 2));
    # Replace the NaN value in the raw dataset (usually -9999) with NaN. Using
    # NaN to keep it consistent with the tsunami data. Also flip the
    # coordinates so it's not upside down (necessary because the asc files
    # store coordinates with an origin of top left, not bottom left)
    z = map(x -> x == missing_val ? NaN : x, reverse(dataset[:, :, 1], dims=2));
    origin, z, x_step, y_step
end

function set_elevations!(network, filepath)
    origin, elevᶻ, x_step, y_step = read_asc(network, filepath);
    max_x_index = size(elevᶻ, 1);
    max_y_index = size(elevᶻ, 2);
    for (id, pos) in network.nodes
        x_index::Int = floor((pos.east - origin[1]) / x_step) + 1;
        y_index::Int = floor((pos.north - origin[2]) / y_step) + 1;
        z = begin
            if 1 ≤ x_index ≤ max_x_index && 1 ≤ y_index ≤ max_y_index
                elevᶻ[x_index,y_index]
            else
                NaN
            end
        end
        network.nodes[id] = OSMX.ENU(pos.east, pos.north, z);
    end
end

function get_slopes(network, segments)
    Dict(map(Iterators.flatten(map(seg -> zip(seg, seg[2:end]), values(segments)))) do ids
        i₁ = network.nodes[ids[1]];
        i₂ = network.nodes[ids[2]];
        rise = i₂.up - i₁.up;
        run = OSMX.distance(i₁, OSMX.ENU(i₂.east, i₂.north, i₁.up));
        slope = rise / run;
        if isnan(slope)
            slope = 0.0;
        end
        (ids, slope)
    end)
end

"""
The mode of a collection; if tied, return the larger value.
"""
function mode(a)
    uniques = reverse(sort(unique(a)));
    counts = map(x -> count(i -> i == x, a), uniques);
    uniques[argmax(counts)]
end

function get_landslides(network, segments, filepath)
    origin, landᶻ, x_step, y_step = read_asc(network, filepath);
    max_x_index = size(landᶻ, 1);
    max_y_index = size(landᶻ, 2);
    coeffs = Dict(
        0 => 1.0,
        1 => 0.5556,
        #2 => 0.5556,
        #3 => 0.5556,
        #4 => 0.5556
    );
    Dict(map(Iterators.flatten(map(seg -> zip(seg, seg[2:end]), values(segments)))) do ids
        p₁ = network.nodes[ids[1]];
        p₂ = network.nodes[ids[2]];
        n = max(round(OSMX.distance(p₁, p₂) / 10m), 1);
        Δp = (p₂ - p₁) / float(n);
        points = map(i -> OSMX.ENU(p₁.east + i*Δp.east, p₁.north + i*Δp.north), 1:n);
        landslides = map(points) do p
            x_index::Int = floor((p.east - origin[1]) / x_step) + 1;
            y_index::Int = floor((p.north - origin[2]) / y_step) + 1;
            if 1 ≤ x_index ≤ max_x_index && 1 ≤ y_index ≤ max_y_index
                landᶻ[x_index,y_index]
            else
                NaN
            end
        end;
        (ids, get(coeffs, mode(landslides), 1.0))
    end)
end

# Get all of the tsunami inundation datasets as a dictionary;
# key is filename besides extension, value is data
datasets = Dict();
geotransform, missing_data_val = nothing, nothing;

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
    tsunami_offset = (AG.getx(network_point, 0), AG.gety(network_point, 0));

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
            map(x -> x[i], AG.applygeotransform.(tuple(geotransform), transform_order...)) .- tsunami_offset[i]
        end
        tsunami_coord(AG.width, 1), reverse(tsunami_coord(AG.height, 2))
    end
    # Replace the NaN value in the raw dataset (usually -9999) with NaN. Can't use
    # `missing` because the `heatmap` function doesn't know how to handle it. Also
    # flip the coordinates so it's not upside down (necessary because the asc files
    # store coordinates with an origin of top left, not bottom left)
    # Use half-mins because that's the interval of the tsunami data
    z = GLMakie.@lift(map(x -> x == missing_data_val ? NaN : x, reverse(datasets[string($half_mins * 30)][:, :, 1], dims=2)));
    x, y, z
end

function tsunami_resolution()::NTuple{2,Float64}
    geotransform[2], -geotransform[6]
end

"""
Find minimum and maximum z values of the tsunami for heatmap normalization.
This can probably be done with `extrema()`, but it's nested so it's
complicated (and it's only computed once, unlike most calculations here)
"""
function tsunami_extrema()::NTuple{2,Float64}
    #minᶻ = minimum(x -> minimum(filter(z -> z != missing_data_val, x[:, :, 1])), values(datasets));
    # Use a min of 0 so ground level is white on the heatmap
    minᶻ = 0.0;
    maxᶻ = maximum(x -> maximum(filter(z -> z != missing_data_val, x[:, :, 1])), values(datasets));
    minᶻ, maxᶻ
end

function people(filename)
    # The people file must have an `Attribute_2` column indicating evacuation or not, regardless of if
    # it will be passed in as a probability when the simulation is run. In that case, it does not matter
    # whether the values in this column are 1 or 0.
    map(CSV.File(filename; normalizenames = true, types = [Int, Float64, Float64, Float64, Bool])) do p
        Person(p.ID, (p.X - offset[1], p.Y - offset[2]), p.Attribute_2)
    end
end

function shelters(network, filename)
    shelter_dict = Dict();
    shelter_locs = CSV.File(filename; normalizenames = true, types = [Int, Float64, Float64]); # Shelter locations
    for shelter in shelter_locs
        node_id = OSMX.nearest_node(network, OSMX.ENU(shelter.x - offset[1], shelter.y - offset[2]));
        shelter_dict[shelter.ID] = Shelter(shelter.ID, node_id, enu_to_tuple(network.nodes[node_id]), false);
    end
    shelter_dict
end

function refresh_data(initial_time, half_mins, network_loc, elev_loc, tsunami_loc, people_loc, shelters_loc)
    new_network, network_segments = network(network_loc);
    set_elevations!(new_network, elev_loc);
    new_network,
    get_slopes(new_network, network_segments),
    tsunami_data(initial_time, half_mins, new_network, tsunami_loc),
    tsunami_extrema(),
    people(people_loc),
    shelters(new_network, shelters_loc)
end

end # module
