```@raw html

<p>The static version of <b>PlanckFunctions.jl</b> Pluto notebook is available at :<a href="../planck_test_git.html">planck_test_git.html</a>
.</p
```


Notebook file with `PlanckFunctions.jl` usage examples are available at [Pluto notebooks](https://github.com/Manarom/PlanckFunctions.jl/blob/master/notebooks).

To run these notebooks, you need:
1) Install `julia` language itself from its official [download page](https://julialang.org/downloads) 
2) Install [Pluto](https://plutojl.org/) notebook from `julia` REPL by entering the following commands line-by-line:
```julia
import Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```
The last line will launch the Pluto starting page in your default browser 
3) Copy/Clone the entire GitHub [repository](https://github.com/Manarom/PlanckFunctions.jl.git) to your local folder. As far as `PlanckFunctions.jl` is not a registered package, all files needed to run the notebooks must be in the `../src` folder with respect to  `note_book_name.jl` file location.
4) Open notebook's `.jl`- file in `Pluto` by providing the full path to the *"Open a notebook"* text field on `Pluto`'s starting page. As far as `Pluto` has its own package manager, it will automatically install all necessary dependancies, which are marked in `using` cell of this file. 
