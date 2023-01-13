pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add CubedSpheres.jl to environment stack

using
  Documenter,
  Literate,
  CairoMakie,
  CubedSpheres

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://clima.github.io/CubedSpheres.jl/stable/",
)

pages = [
    "Home" => "index.md",
    "Conformal Cubed Sphere" => "conformal_cubed_sphere.md",
    "Library" => [ 
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],
]

makedocs(
   sitename = "CubedSpheres.jl",
    modules = [CubedSpheres],
     format = format,
      pages = pages,
    doctest = true,
     strict = true,
      clean = true,
  checkdocs = :exports
)

deploydocs(        repo = "github.com/CliMA/CubedSpheres.jl.git",
               versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
              forcepush = true,
              devbranch = "main",
           push_preview = true
)
