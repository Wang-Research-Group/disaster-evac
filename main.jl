using OpenStreetMapX;
using CSV;
using Makie;
using AbstractPlotting;
using AbstractPlotting.MakieLayout;
using Agents;
using Distributions;
using ArchGDAL;
const OSMX = OpenStreetMapX;
const AP = AbstractPlotting;
const AG = ArchGDAL;


# Unit conversions
const s = 1.0; # Number of frames per second
const ft = 0.3048; # ft to m
const mins = 60s; # mins to frames
const hr = 60mins; # hrs to frames
const mi = 5280ft; # miles to m
const mph = 1mi / 1hr; # mph to m/frame
const m = 1.0 # Exactly 1, no conversion (included for completeness)


road_data = OSMX.get_map_data("data/map.osm", use_cache = false, trim_to_connected_graph = true);

enu_to_tuple(p::OSMX.ENU) = (p.east, p.north);

edges = map(p -> (enu_to_tuple(road_data.nodes[p[1]]), enu_to_tuple(road_data.nodes[p[2]])), road_data.e);

# Set up observables
initial_time = convert(Int, 0s); # Value in frames
curr_time = Makie.Node{Int}(initial_time); # Value in frames
half_mins = Makie.Node{Int}(floor(initial_time / 30s)); # Value rounded down in half-minutes

# Set up stats
num_evacuated = 0;
num_dead = 0;

# Get all of the tsunami inundation datasets as a dictionary;
# key is filename besides extension, value is data
datasets = Dict();
dir_name = "data/tsunami_inundation";
ext = ".asc";
for filename in filter(x -> endswith(x, ext), readdir(dir_name))
    datasets[filename[1:end-length(ext)]] = AG.readraster(dir_name * "/" * filename)
end
datasets[string(initial_time)] = datasets["30"];

dataset = datasets[string(initial_time)];
source = AG.importWKT(AG.getproj(dataset)); # For the tsunami asc files, the projection is UTM
lla = AG.importEPSG(4326); # Lat-long
geotransform = AG.getgeotransform(dataset);

# Determine how far off the tsunami dataset origin is from the road network origin, as ENU/UTM
road_center = OSMX.center(road_data.bounds); # Lat-long of center of the road network, i.e., (0, 0) ENU
network_point = AG.createpoint(road_center.lat, road_center.lon);
AG.createcoordtrans(lla, source) do transform
    AG.transform!(network_point, transform)
end
tsunami_offset = (AG.getx(network_point, 0), AG.gety(network_point, 0));

missing_data_val = AG.getnodatavalue(AG.getband(datasets[string(initial_time)], 1));
# Determine tsunami coordinates in ENU with the same origin as the road
# network. UTM and ENU are both in meters, so no additional conversion needed
tsunamiˣ, tsunamiʸ = begin
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
tsunamiᶻ = @lift(map(x -> x == missing_data_val ? NaN : x, reverse(datasets[string($half_mins * 30)][:, :, 1], dims=2)));
# This can probably be done with `extrema()`, but it's nested so it's
# complicated (and it's only computed once, unlike most calculations here)
minᶻ = minimum(x -> minimum(filter(z -> z != missing_data_val, x[:, :, 1])), values(datasets));
maxᶻ = maximum(x -> maximum(filter(z -> z != missing_data_val, x[:, :, 1])), values(datasets));

"""
The height of the tsunami at a given position
"""
function tsunami_height(position::NTuple{2,Float64})
    x_step = geotransform[2];
    y_step = -geotransform[6];
    x_index::Int = floor((position[1] - tsunamiˣ[1]) / x_step) + 1;
    y_index::Int = floor((position[2] - tsunamiʸ[1]) / y_step) + 1;
    tsunamiᶻ.val[x_index,y_index]
end


people_data = CSV.File("data/pop_coordinates.csv"; normalizenames = true, types = [Int, Float64, Float64, Float64, Bool]);
num_residents = length(people_data); # Number of agents
shelter_locs = CSV.File("data/shelter_coordinates.csv"; normalizenames = true, types = [Int, Float64, Float64]); # Shelter locations
shelters = [begin
        node_id = OSMX.nearest_node(road_data, OSMX.ENU(shelter.x - tsunami_offset[1], shelter.y - tsunami_offset[2]));
        (node_id, enu_to_tuple(road_data.nodes[node_id]))
    end
    for shelter in shelter_locs];

# Time to initialize the agent model

mutable struct Resident <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64} # Vector stored as (x, y) components
    θ::Real # Orientation
    speed::Float64
    dest::Tuple{Int64,NTuple{2,Float64}} # (Intersection ID, (x coord, y coord))
    evacuating::Bool
    alive::Bool
    time_remaining::Int # Number of frames remaining until destination
    milling::Float64 # The amount of time the resident waits until leaving (frames)
end

"""
Initializes a `Resident` if they will never evacuate; several attributes become irrelevant
"""
function Resident(id, pos, θ, alive)
    Resident(id, pos, (0.0, 0.0), θ, 0.0, (0, (0.0, 0.0)), false, alive, typemax(Int), Inf)
end

mutable struct Pedestrian <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64} # Vector stored as (x, y) components
    θ::Real # Orientation
    speed::Float64
    dest::Tuple{Int64,NTuple{2,Float64}} # (Intersection ID, (x coord, y coord))
    path::Array{Tuple{Int64,NTuple{2,Float64}},1}
    alive::Bool
    time_remaining::Int # Number of frames remaining until destination
end

mutable struct Car <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64} # Vector stored as (x, y) components
    θ::Real # Orientation
    speed::Float64
    dest::Tuple{Int64,NTuple{2,Float64}} # (Intersection ID, (x coord, y coord))
    path::Array{Tuple{Int64,NTuple{2,Float64}},1}
    alive::Bool
    ahead::Union{Nothing,Car} # The car ahead
    behind::Union{Nothing,Car} # The car behind (used for optimizations)
    update_ahead::Bool # If the car ahead needs to be updated (used for optimizations)
end

"""
Max that an agent can overshoot or undershoot the destination due to a discrete timestep
"""
discrete_err(a::AbstractAgent) = 0.95 * a.speed;

"""
Velocity vector based on an angle and speed
"""
velocity(θ::Real, speed::Float64) = speed .* (cos(θ), sin(θ));

"""
Update an agent's velocity vector based on its angle and speed
"""
function update_vel!(a::AbstractAgent)
    a.vel = velocity(a.θ, a.speed);
end

"""
Distance between two points
"""
dist(a::NTuple{2,Float64}, b::NTuple{2,Float64}) = distance(OSMX.ENU(a...), OSMX.ENU(b...));

# Determine the angle between two points
angle(a::NTuple{2,Float64}, b::NTuple{2,Float64}) = atan(b[2] - a[2], b[1] - a[1]);

# Angle between the agent and its destination
dest_θ(a::AbstractAgent) = angle(a.pos, a.dest[2]);

# Make the agent face its destination
function face_dest!(a::AbstractAgent)
    a.θ = dest_θ(a);
    update_vel!(a);
end

function get_path(curr_node, safety_node)
    # Only the node IDs
    node_path = OSMX.shortest_route(road_data, curr_node, safety_node)[1][2:end];
    map(x -> (x, enu_to_tuple(road_data.nodes[x])), node_path)
end

next_dest(path) = popfirst!(path);

function next_dest!(a::Union{Pedestrian,Car})
    a.dest = next_dest(a.path);
end

# Color of the agent on the plot
function color(a::AbstractAgent)
    if !a.alive
        return :red;
    end
    if a isa Pedestrian
        return :orange;
    elseif a isa Car
        return :yellow;
    end
    :green
end

# Marker of the agent on the plot
marker(a::AbstractAgent) = a.alive ? :circle : :x;

# Check if the tsunami killed the agent
function tsunami_killed!(agent::AbstractAgent)
    if !agent.alive
        return true;
    end
    if tsunami_height(agent.pos) ≥ 0.5m
        set_speed!(agent, 0.0);
        agent.alive = false;
        global num_dead += 1;
        return true;
    end
    false
end

function agent_step!(resident::Resident, model)
    if tsunami_killed!(resident) || curr_time.val < resident.milling
        return; # Done with the function
    end
    if resident.time_remaining ≤ 0
        # On the off-chance the intersection is safety
        if resident.dest[1] ∈ map(shelter -> shelter[1], shelters)
            Agents.kill_agent!(resident, model);
            global num_evacuated += 1;
            return;
        end
        # Kill the resident, replace with a car or pedestrian
        new_agent_type = rand((Car, Pedestrian));
        new_agent = if new_agent_type == Pedestrian
            new_pedestrian(resident.dest, model)
        else
            new_car(resident.dest, model)
        end
        Agents.kill_agent!(resident, model);
        Agents.add_agent_pos!(new_agent, model);
        if new_agent_type == Car
            # Determine the car ahead; need to do it now, after the agent has been added to the model
            ahead = car_ahead(new_agent, model);
            new_agent.ahead = ahead;
        end
        return;
    end
    Agents.move_agent!(resident, model);
    resident.time_remaining -= 1;
end

function agent_step!(pedestrian::Pedestrian, model)
    if tsunami_killed!(pedestrian)
        return; # Done with the function
    end
    if pedestrian.time_remaining ≤ 0
        pedestrian.pos = pedestrian.dest[2];
        if isempty(pedestrian.path)
            # Reached safety
            Agents.kill_agent!(pedestrian, model);
            global num_evacuated += 1;
            return;
        end
        next_dest!(pedestrian);
        face_dest!(pedestrian);
        update_time_remaining!(pedestrian);
        # If the next intersection is very close, just jump to that instead of a standard move
        if pedestrian.time_remaining ≤ 0
            pedestrian.pos = pedestrian.dest[2];
            return;
        end
    end
    Agents.move_agent!(pedestrian, model);
    pedestrian.time_remaining -= 1;
end

function agent_step!(car::Car, model)
    if tsunami_killed!(car)
        return; # Done with the function
    end
    if dist(car.pos, car.dest[2]) ≤ discrete_err(car)
        # Since speed was already updated, don't have to worry about a vehicle in front
        car.pos = car.dest[2];
        # Tell the car behind that there's nothing ahead anymore
        if !isnothing(car.behind)
            car.behind.ahead = nothing;
        end
        if isempty(car.path)
            # Reached safety
            Agents.kill_agent!(car, model);
            global num_evacuated += 1;
            return;
        end
        next_dest!(car);
        face_dest!(car);
        # During the next decision-making process, update car ahead
        car.update_ahead = true;
        # Since the car jumped to the intersection, don't make a move as well
        return;
    end
    Agents.move_agent!(car, model);
end

"""
Actuates before `agent_step!`, used for setting the car following parameters
"""
function model_step!(model)
    for car ∈ collect(values(filter(p -> p.second isa Car, model.agents)))
        update_speed!(car, model);
    end
end

"""
Generates a random point in the range of the tsunami data; not used currently
"""
function random_point()
    (rand() * (tsunamiˣ[end] - tsunamiˣ[1]) + tsunamiˣ[1],
    rand() * (tsunamiʸ[end] - tsunamiʸ[1]) + tsunamiʸ[1])
end

"""
Free road term of the IDM car following model
"""
function idmᶠʳᵉᵉ(v, v₀)
    a = 0.73m/s^2; # Maximum acceleration
    δ = 4; # Complexity (acceleration exponent)
    a*(1 - (v / v₀)^δ)
end

"""
Interaction term of the IDM car following model
"""
function idmⁱⁿᵗ(v, Δv, sₐ)
    a = 0.73m/s^2; # Maximum acceleration
    b = 1.67m/s^2; # Comfortable braking deceleration (positive)
    s₀ = 2m; # Minimum desired net distance. Can't move if this distance is not met.
    T = 1.5s; # Desired time headway
    sₒₚₜ(v, Δv) = s₀ + v*T + v*Δv / (2*√(a*b)); # Desired gap
    -a*(sₒₚₜ(v, Δv) / sₐ)^2
end

"""
Acceleration produced by the IDM car following model when there is nothing to follow
"""
idm(v, v₀) = idmᶠʳᵉᵉ(v, v₀);

"""
Acceleration produced by the IDM car following model when there is a car to follow
"""
idm(v, v₀, Δv, sₐ) = idmᶠʳᵉᵉ(v, v₀) + idmⁱⁿᵗ(v, Δv, sₐ);

function car_ahead(car, model)
    rotation = [cos(car.θ) sin(car.θ); -sin(car.θ) cos(car.θ)]; # Rotation matrix for clockwise rotation by θ
    # Get all vehicles within radius of next intersection
    max_dist = dist(car.pos, car.dest[2]);
    neighbors = filter(x -> x isa Car, Agents.space_neighbors(car, model, max_dist));
    # Filter for vehicles in front with similar orientation
    nearby = filter(map(neighbors) do id
        (id, rotation * collect(model[id].pos .- car.pos), model[id].θ)
    end) do a
        a[2][1] > 0 && abs(a[2][2]) < 1m && abs(car.θ - a[3]) ≤ π/2
    end
    if isempty(nearby)
        return nothing;
    end
    # Find closest of those vehicles (smallest x value after rotation)
    model[nearby[findmin(map(a -> a[2][1], nearby))[2]][1]]
end

"""
Set the agent's speed to a value and update its velocity
"""
function set_speed!(a::AbstractAgent, speed)
    a.speed = speed;
    update_vel!(a);
end

function update_speed!(car::Car, model)
    # Get the car ahead, up until the next intersection
    if car.update_ahead
        car.ahead = car_ahead(car, model);
        if !isnothing(car.ahead)
            car.ahead.behind = car;
        end
    end
    v₀ = 25mph; # Assume all road speed limits of 25 mph
    l = 20.0ft; # Length of a vehicle
    # Get the acceleration of the vehicle
    a = isnothing(car.ahead) ? idm(car.speed, v₀) : idm(car.speed, v₀, car.speed - car.ahead.speed, dist(car.pos, car.ahead.pos) - l);
    # Update the speed based on the acceleration
    new_speed = max(car.speed + a, 0); # Minimum speed is 0
    set_speed!(car, new_speed);
end

function update_time_remaining!(a::AbstractAgent)
    a.time_remaining = time_remaining(a.pos, a.dest, a.speed);
end

function init_model()
    space = Agents.ContinuousSpace(2; periodic = false); # 2D space

    properties = Dict();
    # The following distributions are determined from a survey
    properties[:ped_shelter_distribution] = Distributions.Gamma(1.920, 1/0.002);
    properties[:car_shelter_distribution] = Distributions.Gamma(1.646, 1/0.000573);

    # Resident needs to go last in each step, otherwise when they transform the agent will get an extra movement
    model = Agents.AgentBasedModel(
        Union{Resident,Pedestrian,Car},
        space;
        scheduler = by_type((Pedestrian, Car, Resident), false),
        properties = properties,
        warn = false
    );

    # Generate a log normal distribution for milling time
    μ = 2.07;
    σ = 0.85;
    milling_distribution = Distributions.LogNormal(μ, σ);

    for person in people_data
        # A distribution provides the number of minutes to mill around (plus the min wait)
        min_wait = 0mins; # Every resident has a milling time of at least this value
        milling_time = rand(milling_distribution)*mins + min_wait;
        resident = new_resident(person, milling_time, model);
        Agents.add_agent_pos!(resident, model);
    end
    Agents.index!(model);
    model
end

"""
Determine number of frames remaining until the agent reaches its destination
"""
time_remaining(pos, dest, speed) = round(Int, dist(pos, dest[2]) / speed);

"""
Like `reduce`, but at each step of accumulation, store it in the list; destructive
"""
function reducemap!(op, list; init = 0)
    list[1] = op(list[1], init);
    for (i, item) in enumerate(list[2:end])
        list[i + 1] = op(item, list[i]);
    end
    list
end

function select_shelter(pos::NTuple{2,Float64}, distribution)
    # Plug in the distance to each shelter to the distribution
    Y = Distributions.pdf.(distribution, dist.(tuple(pos), map(x -> x[2], shelters)));
    # Normalize based on Luce's choice axiom
    shelter_probs = Y ./ sum(Y); # Softmax function; don't need exp because they're already probabilities
    # Sum the probabilities so they go from 0-1
    reducemap!(+, shelter_probs);
    # Determine the index the random number would fall between, then return the corresponding shelter ID
    shelters[searchsortedfirst(shelter_probs, rand())][1]
end

function new_resident(person, milling_time, model)
    id = Agents.nextid(model);
    #id = person.ID;
    #pos = random_point();
    pos = (person.X - tsunami_offset[1], person.Y - tsunami_offset[2]);
    will_evac = person.Attribute_2;
    # Initialize resident facing to the right
    θ = 0;
    resident = if will_evac
        speed = 5mph;
        vel = velocity(θ, speed);
        dest_id = OSMX.nearest_node(road_data, OSMX.ENU(pos...));
        dest = (dest_id, enu_to_tuple(road_data.nodes[dest_id]));
        remaining_time = time_remaining(pos, dest, speed);
        Resident(id, pos, vel, θ, speed, dest, false, true, remaining_time, milling_time)
    else
        Resident(id, pos, θ, true)
    end
    face_dest!(resident);
    resident
end

"""
Initialize a new Pedestrian
"""
function new_pedestrian(resident_pos, model)
    id = Agents.nextid(model);
    speed = 5mph;
    # Initialize pedestrian facing to the right
    θ = 0;
    vel = velocity(θ, speed);
    shelter = select_shelter(resident_pos[2], model.ped_shelter_distribution);
    path = get_path(resident_pos[1], shelter);
    dest = next_dest(path);
    # resident_pos and dest format are (intersection ID, (x coord, y coord))
    pos = resident_pos[2];
    remaining_time = time_remaining(pos, dest, speed);
    pedestrian = Pedestrian(id, pos, vel, θ, speed, dest, path, true, remaining_time);
    # Turn to face the next intersection
    face_dest!(pedestrian);
    pedestrian
end

"""
Initialize a new Car
"""
function new_car(resident_pos, model)
    id = Agents.nextid(model);
    speed = 25mph;
    # Initialize car facing to the right
    θ = 0;
    vel = velocity(θ, speed);
    shelter = select_shelter(resident_pos[2], model.car_shelter_distribution);
    path = get_path(resident_pos[1], shelter);
    dest = next_dest(path);
    # resident_pos and dest format are (intersection ID, (x coord, y coord))
    car = Car(id, resident_pos[2], vel, θ, speed, dest, path, true, nothing, nothing, true);
    # Turn to face the next intersection
    face_dest!(car);
    car
end

model = init_model();

# Lists of tuples `(current time, statistic)`
evac_list = Makie.Node{Array{Tuple{Int,Int},1}}([(0, 0)]);
death_list = Makie.Node{Array{Tuple{Int,Int},1}}([(0, 0)]);

# Get a new agent list (with only enough info for plotting) and step the model every time there's a new frame
agent_list = lift(curr_time; typ = Array{Tuple{Tuple{Float64,Float64},Symbol,Symbol},1}) do now
    # Things to do on every frame
    minute = convert(Int, mins); # Make sure a minute can be evenly divided by frames
    if now % minute == 0
        evac_list[] = push!(evac_list[], (now / minute, num_evacuated));
        death_list[] = push!(death_list[], (now / minute, num_dead));
    end
    curr_half_mins::Int = floor(now / 30s);
    # Need to do a conditional because we don't want the observable to trigger unless it's a new value
    if curr_half_mins ≠ half_mins.val
        half_mins[] = curr_half_mins;
    end
    if now ≠ initial_time
        # Step, with the model step happening before the agent step
        Agents.step!(model, agent_step!, model_step!, 1, false);
    end
    map(a -> (a.pos, color(a), marker(a)), Agents.allagents(model))
end

# Separate out the tuple of info into new observables because `scatter` can't handle it otherwise
position_list = @lift(getindex.($agent_list, 1));
color_list = @lift(getindex.($agent_list, 2));
marker_list = @lift(getindex.($agent_list, 3));

#=
hour = convert(Int, hr);
for now in 1:hour
    curr_time[] = now
    sleep(.0001)
end
println(num_evacuated);
println(num_dead);
=#

scene, layout = MakieLayout.layoutscene(resolution = (1200, 900))
main_scene = layout[1:2, 1] = MakieLayout.LAxis(scene)
evac_plot = layout[1, 2] = MakieLayout.LAxis(scene, xlabel = "Minutes", ylabel = "Total Evacuated", title = "Successful Evacuations")
death_plot = layout[2, 2] = MakieLayout.LAxis(scene, xlabel = "Minutes", ylabel = "Total Deaths", title = "Deaths")
button = layout[0, :] = MakieLayout.LButton(scene, label = "Start/Stop")
button.tellwidth = false

AbstractPlotting.linesegments!(main_scene, edges)

# Render a heatmap
AbstractPlotting.heatmap!(main_scene, tsunamiˣ, tsunamiʸ, tsunamiᶻ; colormap = :GnBu_9, colorrange = (minᶻ, maxᶻ))

AbstractPlotting.scatter!(main_scene, position_list; color = color_list, markersize = 5, marker = marker_list);

AbstractPlotting.scatter!(evac_plot, evac_list);
MakieLayout.limits!(evac_plot, 0, 60, 0, num_residents);

AbstractPlotting.scatter!(death_plot, death_list);
MakieLayout.limits!(death_plot, 0, 60, 0, num_residents);

AbstractPlotting.scatter!(main_scene, [shelter[2] for shelter in shelters]; color = :blue);

on(button.clicks) do click
    if click == 1
        hour = convert(Int, hr); # Make sure an hour can be evenly divided by frames
        @async for now in 1:hour
            while button.clicks.val % 2 == 0
                sleep(.5)
            end
            curr_time[] = now
            sleep(.0001)
            if now == hour
                button.clicks[] = 0;
                println(num_evacuated);
                println(num_dead);
                global model = init_model();
                curr_time[] = initial_time;
                evac_list[] = [(0, 0)];
                global num_evacuated = 0;
                death_list[] = [(0, 0)];
                global num_dead = 0;
            end
        end
    end
end

main_scene.aspect = MakieLayout.DataAspect()
MakieLayout.hidedecorations!(main_scene)

scene