module DisasterEvac

import OpenStreetMapX;
import CSV;
import Makie;
import Agents;
import Distributions;
import ArchGDAL;
const OSMX = OpenStreetMapX;
const AG = ArchGDAL;

include("utilities.jl");
using .Utils;

include("parse_data.jl");


const Agent = Agents.AbstractAgent;

# Speed and min milling time
const default_params = (25mph, 0mins);
const death_threshold = 0.5m;


# Set up observables
initial_time = convert(Int, 0s); # Value in frames
curr_time = Makie.Node{Int}(initial_time); # Value in frames
half_mins = Makie.Node{Int}(floor(initial_time / 30s)); # Value rounded down in half-minutes

# Set up stats
evacuated = [];
dead = [];

function set_filenames()
    filenames = Dict();
    filenames["network"] = "src/data/map.osm";
    filenames["elevation"] = "src/data/elevation/elevation.asc";
    filenames["tsunami"] = "src/data/tsunami_inundation";
    filenames["people"] = "src/data/pop_coordinates.csv";
    filenames["shelters"] = "src/data/shelter_coordinates.csv";
    filenames
end

filenames = set_filenames();


network, network_segments = Data.network(filenames["network"]);

Data.set_elevations!(network, filenames["elevation"]);
slopes = Data.get_slopes(network, network_segments);

tsunamiˣ, tsunamiʸ, tsunamiᶻ = Data.tsunami_data(initial_time, half_mins, network, filenames["tsunami"]);
minᶻ, maxᶻ = Data.tsunami_extrema();

people = Data.people(filenames["people"]);
num_residents = length(people); # Number of agents
shelters = Data.shelters(network, filenames["shelters"]);

function refresh_data()::Nothing
    data = Data.refresh_data(initial_time, half_mins,
                             filenames["network"], filenames["elevation"], filenames["tsunami"],
                             filenames["people"], filenames["shelters"]
    );
    network = data[1];
    slopes = data[2];
    tsunamiˣ, tsunamiʸ, tsunamiᶻ = data[3];
    minᶻ, maxᶻ = data[4];
    people = data[5];
    shelters = data[6];
    nothing
end


# Time to initialize the agent model

mutable struct Resident <: Agent
    id::Int
    pos::Point
    vel::Point # Vector stored as (x, y) components
    θ::Real # Orientation
    speed::Float64
    dest::Tuple{Int64,Point} # (Intersection ID, (x coord, y coord))
    evacuating::Bool
    alive::Bool
    time_remaining::Int # Number of frames remaining until destination
    milling::Float64 # The amount of time the resident waits until leaving (frames)
    ext_id::Int # External ID, used for results
end

"""
Initializes a `Resident` if they will never evacuate; several attributes become irrelevant.
"""
function Resident(id, pos, θ, alive, ext_id)::Resident
    Resident(id, pos, (0.0, 0.0), θ, 0.0, (0, (0.0, 0.0)), false, alive, typemax(Int), Inf, ext_id)
end

mutable struct Pedestrian <: Agent
    id::Int
    pos::Point
    vel::Point # Vector stored as (x, y) components
    θ::Real # Orientation
    speed::Float64
    dest::Tuple{Int64,Point} # (Intersection ID, (x coord, y coord))
    path::Array{Tuple{Int64,Point},1}
    alive::Bool
    time_remaining::Int # Number of frames remaining until destination
    ext_id::Int # External ID, used for results
end

mutable struct Car <: Agent
    id::Int
    pos::Point
    vel::Point # Vector stored as (x, y) components
    θ::Real # Orientation
    speed::Float64
    dest::Tuple{Int64,Point} # (Intersection ID, (x coord, y coord))
    path::Array{Tuple{Int64,Point},1}
    alive::Bool
    ahead::Union{Nothing,Car} # The car ahead
    behind::Union{Nothing,Car} # The car behind (used for optimizations)
    update_ahead::Bool # If the car ahead needs to be updated (used for optimizations)
    ext_id::Int # External ID, used for results
end

"""
Max that an agent can overshoot or undershoot the destination due to a discrete timestep.
"""
discrete_err(a::Agent)::Float64 = 0.95 * a.speed;

"""
Velocity vector based on an angle and speed.
"""
velocity(θ::Real, speed::Float64)::Point = speed .* (cos(θ), sin(θ));

"""
Update an agent's velocity vector based on its angle and speed.
"""
function update_vel!(a::Agent)::Nothing
    a.vel = velocity(a.θ, a.speed);
    nothing
end

"""
Angle between the agent and its destination.
"""
dest_θ(a::Agent)::Float64 = Utils.angle(a.pos, a.dest[2]);

"""
Make the agent face its destination.
"""
function face_dest!(a::Agent)::Nothing
    a.θ = dest_θ(a);
    update_vel!(a);
end

"""
Determine the shortest path to the shelter.
"""
function get_path(curr_node, shelter, segments)::Array{Tuple{Int64,Point},1}
    # Only the node IDs
    node_path = OSMX.shortest_route(network, curr_node, shelters[shelter].node_id)[1];
    map(x -> (x, enu_to_tuple(network.nodes[x])), Iterators.flatten(map(i -> segments[i][2:end], zip(node_path, node_path[2:end]))))
end

"""
Find the next step of the path.
"""
next_dest(path)::Tuple{Int64,Point} = popfirst!(path);

"""
Advance to the next step of the path.
"""
function next_dest!(a::Union{Pedestrian,Car})::Nothing
    a.dest = next_dest(a.path);
    nothing
end

"""
Color of the agent on the plot.
"""
function color(a::Agent)::Symbol
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

"""
Marker of the agent on the plot.
"""
marker(a::Agent)::Symbol = a.alive ? :circle : :x;

"""
The height of the tsunami at a given position.
"""
function tsunami_height(position::Point)::Float64
    x_step, y_step = Data.tsunami_resolution();
    x_index::Int = floor((position[1] - tsunamiˣ[1]) / x_step) + 1;
    y_index::Int = floor((position[2] - tsunamiʸ[1]) / y_step) + 1;
    tsunamiᶻ.val[x_index,y_index]
end

"""
Check if the tsunami killed the agent.
"""
function tsunami_killed!(agent::Agent)::Bool
    if !agent.alive
        return true;
    end
    if tsunami_height(agent.pos) ≥ death_threshold
        set_speed!(agent, 0.0);
        agent.alive = false;
        push!(dead, (agent.ext_id, agent.pos));
        return true;
    end
    false
end

"""
Advance the step of a resident.
"""
function agent_step!(resident::Resident, model)::Nothing
    if tsunami_killed!(resident) || curr_time.val < resident.milling
        return; # Done with the function
    end
    if resident.time_remaining ≤ 0
        # On the off-chance the intersection is safety
        shelter_id = findfirst(shelter -> shelter.node_id == resident.dest[1], shelters);
        if !isnothing(shelter_id)
            push!(evacuated, (resident.ext_id, shelter_id));
            Agents.kill_agent!(resident, model);
            return;
        end
        # Kill the resident, replace with a car or pedestrian
        new_agent_type = rand((Car, Pedestrian));
        new_agent = if new_agent_type == Pedestrian
            new_pedestrian(resident.ext_id, resident.dest, model)
        else
            new_car(resident.ext_id, resident.dest, model)
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
    nothing
end

"""
Advance the step of a pedestrian.
"""
function agent_step!(pedestrian::Pedestrian, model)::Nothing
    if tsunami_killed!(pedestrian)
        return; # Done with the function
    end
    if pedestrian.time_remaining ≤ 0
        pedestrian.pos = pedestrian.dest[2];
        if isempty(pedestrian.path)
            # Reached safety
            push!(evacuated, (pedestrian.ext_id, shelter_at_node(pedestrian.dest[1])));
            Agents.kill_agent!(pedestrian, model);
            return;
        end
        # Update speed
        road_slope = slopes[(pedestrian.dest[1], pedestrian.path[1][1])];
        set_speed!(pedestrian, ped_speed(road_slope));

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
    nothing
end

"""
Advance the step of a car.
"""
function agent_step!(car::Car, model)::Nothing
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
            push!(evacuated, (car.ext_id, shelter_at_node(car.dest[1])));
            Agents.kill_agent!(car, model);
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
    nothing
end

"""
Actuates before `agent_step!`; used for checking for flooded shelters and
setting the car following parameters.
"""
function model_step!(model)::Nothing
    for (id, shelter) ∈ shelters
        if tsunami_height(shelter.pos) ≥ death_threshold
            shelter.inundated = true;
            had_evacuated = findall(x -> x[2] == id, evacuated);
            # Remove dead from evacuated
            now_dead = popat!.(tuple(evacuated), reverse(had_evacuated)); # Have to reverse so popping uses correct indices
            # Add to dead
            push!(dead, map(x -> (x[1], shelter.pos), now_dead)...);
        else
            shelter.inundated = false;
        end
    end

    for car ∈ collect(values(filter(p -> p.second isa Car, model.agents)))
        update_speed!(car, model);
    end
    nothing
end

"""
Generates a random point in the range of the tsunami data; not used currently.
"""
function random_point()::Point
    (rand() * (tsunamiˣ[end] - tsunamiˣ[1]) + tsunamiˣ[1],
    rand() * (tsunamiʸ[end] - tsunamiʸ[1]) + tsunamiʸ[1])
end

ped_speed(slope)::Float64 = 1.65*exp(-2.30*abs(slope - 0.004))*m/s;

"""
Free road term of the IDM car following model.
"""
function idmᶠʳᵉᵉ(v, v₀)::Float64
    a = 0.73m/s^2; # Maximum acceleration
    δ = 4; # Complexity (acceleration exponent)
    a*(1 - (v / v₀)^δ)
end

"""
Interaction term of the IDM car following model.
"""
function idmⁱⁿᵗ(v, Δv, sₐ)::Float64
    a = 0.73m/s^2; # Maximum acceleration
    b = 1.67m/s^2; # Comfortable braking deceleration (positive)
    s₀ = 2m; # Minimum desired net distance. Can't move if this distance is not met.
    T = 1.5s; # Desired time headway
    sₒₚₜ(v, Δv) = s₀ + v*T + v*Δv / (2*√(a*b)); # Desired gap
    -a*(sₒₚₜ(v, Δv) / sₐ)^2
end

"""
Acceleration produced by the IDM car following model when there is nothing to follow.
"""
idm(v, v₀)::Float64 = idmᶠʳᵉᵉ(v, v₀);

"""
Acceleration produced by the IDM car following model when there is a car to follow.
"""
idm(v, v₀, Δv, sₐ)::Float64 = idmᶠʳᵉᵉ(v, v₀) + idmⁱⁿᵗ(v, Δv, sₐ);

"""
Determine the car ahead, up to the next intersection.
"""
function car_ahead(car, model)::Union{Car,Nothing}
    rotation = [cos(car.θ) sin(car.θ); -sin(car.θ) cos(car.θ)]; # Rotation matrix for clockwise rotation by θ
    # Get all vehicles within radius of next intersection
    max_dist = dist(car.pos, car.dest[2]);
    neighbors = filter(x -> model[x] isa Car, Agents.space_neighbors(car, model, max_dist));
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
Set the agent's speed to a value and update its velocity.
"""
function set_speed!(a::Agent, speed)::Nothing
    a.speed = speed;
    update_vel!(a);
    nothing
end

"""
Update the speed of a car.
"""
function update_speed!(car::Car, model)::Nothing
    # Get the car ahead, up until the next intersection
    if car.update_ahead || car.speed == 0
        car.ahead = car_ahead(car, model);
        if !isnothing(car.ahead)
            car.ahead.behind = car;
        end
    end
    car.update_ahead = false;
    v₀ = model.speed_limit;
    l = 20.0ft; # Length of a vehicle
    # Get the acceleration of the vehicle
    a = isnothing(car.ahead) ? idm(car.speed, v₀) : idm(car.speed, v₀, car.speed - car.ahead.speed, dist(car.pos, car.ahead.pos) - l);
    # Update the speed based on the acceleration
    new_speed = max(car.speed + a, 0); # Minimum speed is 0
    set_speed!(car, new_speed);
    nothing
end

"""
Determine number of frames remaining until the agent reaches its next intersection.
"""
time_remaining(pos, dest, speed)::Int64 = round(Int, dist(pos, dest[2]) / speed);

"""
Update the time remaining of an agent.

For more information, refer to `time_remaining()`.
"""
function update_time_remaining!(a::Agent)::Nothing
    a.time_remaining = time_remaining(a.pos, a.dest, a.speed);
    nothing
end

function shelter_at_node(n)
    findfirst(x -> x.node_id == n, shelters)
end

"""
Initialize the model.
"""
function init_model(speed_limit, min_wait)::Agents.ABM
    space = Agents.ContinuousSpace(2; periodic = false); # 2D space

    properties = Dict();
    # The following distributions are determined from a survey
    properties[:ped_shelter_distribution] = Distributions.Gamma(1.920, 1/0.002);
    properties[:car_shelter_distribution] = Distributions.Gamma(1.646, 1/0.000573);

    # Set all road speed limits in mph
    properties[:speed_limit] = speed_limit;

    # Resident needs to go last in each step, otherwise when they transform the agent will get an extra movement
    model = Agents.AgentBasedModel(
        Union{Resident,Pedestrian,Car},
        space;
        scheduler = Agents.by_type((Pedestrian, Car, Resident), false),
        properties = properties,
        warn = false
    );

    # Generate a log normal distribution for milling time
    μ = 2.07;
    σ = 0.85;
    milling_distribution = Distributions.LogNormal(μ, σ);

    for person in people
        # A distribution provides the number of minutes to mill around (plus the min wait)
        # Min wait is in mins
        milling_time = rand(milling_distribution)*mins + min_wait;
        resident = new_resident(person, milling_time, model);
        Agents.add_agent_pos!(resident, model);
    end
    Agents.index!(model);
    model
end

init_model() = init_model(default_params...);

"""
Selects a shelter based on how far away it is and a provided distribution.
"""
function select_shelter(pos::Point, distribution)::Int64
    shelter_ids = collect(keys(shelters));
    # Plug in the distance to each shelter to the distribution
    Y = Distributions.pdf.(distribution, dist.(tuple(pos), map(x -> x.pos, values(shelters))));
    # Normalize based on Luce's choice axiom
    shelter_probs = Y ./ sum(Y); # Softmax function; don't need exp because they're already probabilities
    # Sum the probabilities so they go from 0-1
    cumsum!(shelter_probs, shelter_probs);
    # Determine the index the random number would fall between, then return the corresponding shelter ID
    shelter_ids[searchsortedfirst(shelter_probs, rand())]
end

"""
Initialize a new resident.
"""
function new_resident(person, milling_time, model)::Resident
    id = Agents.nextid(model);
    ext_id = person.id;
    #pos = random_point();
    pos = person.pos;
    #will_evac = person.Attribute_2;
    will_evac = person.evac;
    # Initialize resident facing to the right
    θ = 0;
    resident = if will_evac
        speed = 5mph;
        vel = velocity(θ, speed);
        dest_id = OSMX.nearest_node(network, OSMX.ENU(pos...));
        dest = (dest_id, enu_to_tuple(network.nodes[dest_id]));
        remaining_time = time_remaining(pos, dest, speed);
        Resident(id, pos, vel, θ, speed, dest, false, true, remaining_time, milling_time, ext_id)
    else
        Resident(id, pos, θ, true, ext_id)
    end
    face_dest!(resident);
    resident
end

"""
Initialize a new pedestrian.
"""
function new_pedestrian(ext_id, resident_pos, model)::Pedestrian
    id = Agents.nextid(model);
    shelter = select_shelter(resident_pos[2], model.ped_shelter_distribution);
    path = get_path(resident_pos[1], shelter, network_segments);
    dest = next_dest(path);

    road_slope = slopes[(resident_pos[1], dest[1])];
    speed = ped_speed(road_slope);

    # Initialize pedestrian facing to the right
    θ = 0;
    vel = velocity(θ, speed);
    # resident_pos and dest format are (intersection ID, (x coord, y coord))
    pos = resident_pos[2];
    remaining_time = time_remaining(pos, dest, speed);
    pedestrian = Pedestrian(id, pos, vel, θ, speed, dest, path, true, remaining_time, ext_id);
    # Turn to face the next intersection
    face_dest!(pedestrian);
    pedestrian
end

"""
Initialize a new car.
"""
function new_car(ext_id, resident_pos, model)::Car
    id = Agents.nextid(model);
    speed = model.speed_limit;
    # Initialize car facing to the right
    θ = 0;
    vel = velocity(θ, speed);
    shelter = select_shelter(resident_pos[2], model.car_shelter_distribution);
    path = get_path(resident_pos[1], shelter, network_segments);
    dest = next_dest(path);
    # resident_pos and dest format are (intersection ID, (x coord, y coord))
    car = Car(id, resident_pos[2], vel, θ, speed, dest, path, true, nothing, nothing, true, ext_id);
    # Turn to face the next intersection
    face_dest!(car);
    car
end

model = init_model();

# Lists of tuples `(current time, statistic)`
evac_list = Makie.Node{Array{Tuple{Int,Int},1}}([(0, 0)]);
death_list = Makie.Node{Array{Tuple{Int,Int},1}}([(0, 0)]);

# Get a new agent list (with only enough info for plotting) and step the model every time there's a new frame
agent_list = Makie.lift(curr_time; typ = Array{Tuple{Tuple{Float64,Float64},Symbol,Symbol},1}) do now
    # Things to do on every frame
    minute = convert(Int, mins); # Make sure a minute can be evenly divided by frames
    if now % minute == 0
        evac_list[] = push!(evac_list[], (now / minute, length(evacuated)));
        death_list[] = push!(death_list[], (now / minute, length(dead)));
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
position_list = Makie.@lift(getindex.($agent_list, 1));
color_list = Makie.@lift(getindex.($agent_list, 2));
marker_list = Makie.@lift(getindex.($agent_list, 3));

fig = Makie.Figure(; resolution = (1200, 900));
simulation = fig[1:2, 1] = Makie.Axis(fig);
evac_plot = fig[1, 2] = Makie.Axis(fig; xlabel = "Minutes", ylabel = "Total Evacuated", title = "Successful Evacuations");
death_plot = fig[2, 2] = Makie.Axis(fig; xlabel = "Minutes", ylabel = "Total Deaths", title = "Deaths");
Makie.linkaxes!(evac_plot, death_plot);
button = fig[0, :] = Makie.Button(fig; label = "Start/Stop");
button.tellwidth = false;

lines = collect(Iterators.flatten(map(seg -> vcat(map(x -> enu_to_tuple(network.nodes[x]), seg), [(NaN, NaN)]), values(network_segments))));
Makie.lines!(simulation, lines);

# Render a heatmap
Makie.heatmap!(simulation, tsunamiˣ, tsunamiʸ, tsunamiᶻ; colormap = :GnBu_9, colorrange = (minᶻ, maxᶻ));

Makie.scatter!(simulation, position_list; color = color_list, markersize = 5, marker = marker_list);

Makie.scatter!(evac_plot, evac_list);
Makie.limits!(evac_plot, 0, 60, 0, num_residents);

Makie.scatter!(death_plot, death_list);
Makie.limits!(death_plot, 0, 60, 0, num_residents);

Makie.scatter!(simulation, [shelter.pos for shelter in values(shelters)]; color = :blue);

simulation.aspect = Makie.DataAspect();
Makie.hidedecorations!(simulation);

function reset_model!(speed_limit, min_wait)::Nothing
    global model = init_model(speed_limit, min_wait);
    curr_time[] = initial_time;
    evac_list[] = [(0, 0)];
    global evacuated = [];
    death_list[] = [(0, 0)];
    global dead = [];
    nothing
end

function run_no_gui(times, options)
    # Interpreting options as mph and mins
    hour = convert(Int, hr);
    stats = [];
    for option ∈ options
        println("Option ", option);
        for i in 1:times
            println("Run ", i);
            reset_model!(float(option[1])*mph, float(option[2])*mins);
            for now in 1:hour
                curr_time[] = now;
            end
            println("Evacuated: ", length(evacuated));
            println("Dead: ", length(dead));
            push!(stats, (evacuated, dead));
        end
    end
    reset_model!(default_params...);
    stats
end

run_no_gui(times) = run_no_gui(times, [(default_params[1]/mph, default_params[2]/mins)]);

run_no_gui() = run_no_gui(1);

function run_gui()::Nothing
    Makie.display(fig);
    nothing
end

Makie.on(button.clicks) do click
    if click == 1
        hour = convert(Int, hr); # Make sure an hour can be evenly divided by frames
        @async for now in 1:hour
            while button.clicks.val % 2 == 0
                sleep(.5);
            end
            curr_time[] = now;
            sleep(.0001); # Needed so the frame can render
            if now == hour
                println("Evacuated: ", length(evacuated));
                println("Dead: ", length(dead));
                button.clicks[] = 0;
            end
        end
    end
end

Makie.on(fig.scene.events.window_open) do status
    # If window just closed, reset for a new run
    if !fig.scene.events.window_open.val
        reset_model!(default_params...);
    end
end

function run_record(filename)::Nothing
    hour = convert(Int, hr);
    Makie.record(fig, filename, 1:hour; framerate = 60) do t
        curr_time[] = t;
    end
    println("Evacuated: ", length(evacuated));
    println("Dead: ", length(dead));
    reset_model!(default_params...);
    nothing
end

end # module
