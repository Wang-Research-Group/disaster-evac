# disaster-evac
Simulation that models evacuation from natural disasters.
Currently designed to be tsunami-specific.

# Installation
To run, open a Julia REPL in the directory of this repo, then execute the following steps:

1. `]activate .`
2. `]instantiate`
3. `include("src/DisasterEvac.jl")`
4. Run the simulation
- With a graphical interface: `DisasterEvac.run_gui()`
- Without a graphical interface: `DisasterEvac.run_no_gui()`
- Record an animation: `DisasterEvac.run_record()`

# Suggested Workflow
We suggest you use Visual Studio Code with the Julia extension when navigating the codebase, as it provides some nice Julia-specific suggestions and allows you to run a REPL at the same time as looking at the code.

To do so:
1. Install Julia, such that you can execute it from a command line.
2. Install Visual Studio Code.
3. Add the `Julia` extension (can add extensions from the left sidebar).
4. Navigate to the codebase with File -> Open Folder....
5. Open a terminal if there isn't one already with Terminal -> New Terminal.
6. Execute Julia (e.g., `julia.exe` on a Windows machine). This should open up a Julia REPL.
7. Follow the three steps in the Installation section above.

# Adjusting Locations
To use a different dataset for a new location, replace the files in the `data` folder.
All replaced files need to be the same name as the original, including extension.
The current codebase assumes the tsunami inundation data is in 30-sec increments, for a total of 1 hour.

# Data Sources
This simulation would not have been possible without data and work that has been collected by others. In this section we identify sources of this information.

## Agents.jl \(dependency\)
George Datseris, Ali R. Vahdati, and Timothy C. DuBois. *Agents.jl: A performant and feature-full agent based modelling software of minimal code complexity*. 2021. arXiv: [2101.10072 \[cs.MA\]](https://arxiv.org/abs/2101.10072).

## Elevation data
OGEO

Oregon Geospatial Enterprise Office, 2017. Oregon 10m Digital Elevation Model (DEM).

## Landslide data
DOGAMI

Gabel, L. L. S., O’Brien, F. E., Bauer, J. M., Allan, J. C., 2019. TSUNAMI EVACUATION ANALYSIS OF COMMUNITIES SURROUNDING THE COOS BAY ESTUARY: BUILDING COMMUNITY RESILIENCE ON THE OREGON COAST. Tech. Rep. O-19-07, Oregon Department of Geology and Mineral Industries.

## Tsunami data
DOGAMI

Witter, R. C., Zhang, Y., Wang, K., Priest, G. R., Goldfinger, C., Stimely, L. L., English,J. T., Ferro, P. A., 2011. Simulating Tsunami Inundation at Bandon, Coos County, Oregon, Using Hypothetical Cascadia and Alaska Earthquake Scenarios. Tech. Rep. Special Paper 43, Oregon Department of Geology and Mineral Industries.

## Population data
2019 census estimate, with modifications based on road density.

## Shelter data
Chen Chen, Michael K. Lindell, Haizhong Wang, 2021. “Tsunami Preparedness and Resilience in the Cascadia Subduction Zone: A Multistage Model of Expected Evacuation Decisions and Mode Choice,” International Journal of Disaster Risk Reduction. Vol. 59 (2021), 102244.