# User guide

## Introduction
This package provides constrained Delaunay-Triangulations. It implements an edge-flipping algorithm to restore the Delaunay property and to remove intersecting edges. The user can define constraints as closed polygons and the package automatically marks faces interior to constrained regions.  The package provides various convenicenc methods, such as exporting triangulations as FEM-meshes or computing the centroid of faces.

## Usage
The first step is to create an empty mesh and set an appropriate bound box. `DelaunayMeshes` internally scales all vertex coordinates to fit into an enclosing virtual triangle. The bounding box establishes this scaling. Note that all vertices that are pushed into the triangulation later on must fit into this bounding box.

```
    import DelaunayMeshes

    # create empty mesh and set bounding box
    local mesh = DelaunayMeshes.Mesh()
    DelaunayMeshes.setboundingbox(mesh, [-15.0, 15.0, -15.0, 15.0])
```

Vertices are added to the triangulation by pushing vectors of coordinates.

```
    local seed = rand(UInt32)
    local nvertices = 200
    local points = rand(Float64, nvertices, 2)*30. - 15.
    push!(mesh, points)  
```

We can plot the triangulation using ``

```
    xc, yc = DelaunayMeshes.getdelaunaycoordinates(mesh.tesselation)
    p = Winston.FramedPlot(aspect_ratio=1)
    Winston.add(p, Winston.Curve(xc, yc))
    Winston.add(p, Winston.Points(
      vorVert[mesh.faceLocation,1], vorVert[mesh.faceLocation,2],
      kind="circle", color="red"))
    Winston.savefig(p, "RandomTriangulation.svg")
```
