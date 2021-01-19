# disaster-evac
Simulation that models evacuation from natural disasters.
Currently designed to be tsunami-specific.

# Installation
To run, make this repo your current working directory, then run `include("main.jl")` in a Julia REPL.

# Suggested Workflow
We suggest you use Visual Studio Code with the Julia extension when navigating the codebase, as it provides some nice Julia-specific suggestions and allows you to run a REPL at the same time as looking at the code.

To do so:
1. Install Julia, such that you can execute it from a command line.
2. Install Visual Studio Code.
3. Add the `Julia` extension (can add extensions from the left sidebar).
4. Navigate to the codebase with File -> Open Folder....
5. Open a terminal if there isn't one already with Terminal -> New Terminal.
6. Execute Julia (e.g., `julia.exe` on a Windows machine). This should open up a Julia REPL.
7. Type `include("main.jl")`.

# Adjusting Locations
To use a different dataset for a new location, replace the files in the `data` folder.
All replaced files need to be the same name as the original, including extension.
The current codebase assumes the tsunami inundation data is in 30-sec increments, for a total of 1 hour.
