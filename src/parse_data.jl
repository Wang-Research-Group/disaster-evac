module Data

#export network, people, shelters, tsunami;
export network;

using OpenStreetMapX;
const OSMX = OpenStreetMapX;

network = OSMX.get_map_data("src/data/map.osm", use_cache = false, trim_to_connected_graph = true);

end