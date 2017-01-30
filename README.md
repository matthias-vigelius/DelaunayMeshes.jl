# DelaunayMeshes
This package provides constrained Delaunay-Triangulations. It implements an edge-flipping algorithm to restore the Delaunay property and to remove intersecting edges. The user can define constraints as closed polygons and the package automatically marks faces interior to constrained regions.  The package provides various convenience methods, such as exporting triangulations as FEM-meshes or computing the centroid of faces.

Please consult the [documentation](https://matthias-vigelius.github.io/DelaunayMeshes.jl/latest/UserGuide.html) for details.

[![Build Status](https://travis-ci.org/matthias-vigelius/DelaunayMeshes.jl.svg?branch=master)](https://travis-ci.org/matthias-vigelius/DelaunayMeshes.jl)

[![Coverage Status](https://coveralls.io/repos/matthias-vigelius/DelaunayMeshes.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/matthias-vigelius/DelaunayMeshes.jl?branch=master)

[![codecov.io](http://codecov.io/github/matthias-vigelius/DelaunayMeshes.jl/coverage.svg?branch=master)](http://codecov.io/github/matthias-vigelius/DelaunayMeshes.jl?branch=master)
