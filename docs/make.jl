using Documenter, Flows

makedocs(sitename="Flows.jl",
         modules = [Flows],
         authors="Davide Lasagna",
         pages = ["Home"                  => "index.md",
                  "Quick Start"           => "quickstart.md",
                  "Monitors"              => "monitors.md",
                  "Coupled systems"       => "coupled.md",
                  "Quadrature equations"  => "quadrature.md",
                  "Available Methods"     => "available-methods.md",
                  "Full API"              => "api.md",
                  ])