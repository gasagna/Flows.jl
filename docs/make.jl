using Documenter, Flows

makedocs(sitename="Flows.jl",
         modules = [Flows],
         authors="Davide Lasagna",
         pages = ["Home"                  => "index.md",
                  "Quick Start"           => "quickstart.md",
                  "Monitors"              => "monitors.md",
                  "Coupled systems"       => "coupled.md",
                  "Quadrature equations"  => "quadrature.md",
                  "Solution storages"     => "storage.md",
                  "Advanced features"     => "advanced.md",
                  "Examples"              => "examples.md",
                  "Available methods"     => "available-methods.md",
                  "Full API"              => "api.md",
                  ])

# deploydocs(
    # repo = "github.com/gasagna/Flows.jl.git",
# )