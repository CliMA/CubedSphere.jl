pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add CubedSphere.jl to environment stack

using
  Documenter,
  Literate,
  CairoMakie,
  CubedSphere

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://clima.github.io/CubedSphere.jl/stable/",
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
   sitename = "CubedSphere.jl",
    modules = [CubedSphere],
     format = format,
      pages = pages,
    doctest = true,
     strict = true,
      clean = true,
  checkdocs = :exports
)

deploydocs(        repo = "github.com/CliMA/CubedSphere.jl.git",
               versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
              forcepush = true,
              devbranch = "main",
           push_preview = true
)
