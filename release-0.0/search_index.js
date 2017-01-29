var documenterSearchIndex = {"docs": [

{
    "location": "UserGuide.html#",
    "page": "User Guide",
    "title": "User Guide",
    "category": "page",
    "text": ""
},

{
    "location": "UserGuide.html#User-guide-1",
    "page": "User Guide",
    "title": "User guide",
    "category": "section",
    "text": ""
},

{
    "location": "UserGuide.html#Introduction-1",
    "page": "User Guide",
    "title": "Introduction",
    "category": "section",
    "text": "This package provides constrained Delaunay-Triangulations. It implements an edge-flipping algorithm to restore the Delaunay property and to remove intersecting edges. The user can define constraints as closed polygons and the package automatically marks faces interior to constrained regions.  The package provides various convenience methods, such as exporting triangulations as FEM-meshes or computing the centroid of faces."
},

{
    "location": "UserGuide.html#Unconstrained-triangulations-1",
    "page": "User Guide",
    "title": "Unconstrained triangulations",
    "category": "section",
    "text": "The first step is to create an empty mesh and set an appropriate bound box. DelaunayMeshes internally scales all vertex coordinates to fit into an enclosing virtual triangle. The bounding box establishes this scaling. Note that all vertices that are pushed into the triangulation later on must fit into this bounding box.import DelaunayMeshes\n\n# create empty mesh and set bounding box\nmesh = DelaunayMeshes.Mesh()\nDelaunayMeshes.setboundingbox(mesh, [-15.0, 15.0, -15.0, 15.0])Vertices are added to the triangulation by pushing vectors of coordinates.seed = rand(UInt32)\nnvertices = 200\npoints = rand(Float64, nvertices, 2)*30. - 15.\npush!(mesh, points)  We can plot the triangulation using getdelaunaycoordinatesimport Winston\nxc, yc = DelaunayMeshes.getdelaunaycoordinates(mesh.tesselation)\np = Winston.FramedPlot(aspect_ratio=1)\nWinston.add(p, Winston.Curve(xc, yc))(Image: Triangulation)Note that, internally, DelaunayMeshes scales down all coordinates into a numerically favorable range. The method unscalepoints allows to recover the original coordinates:unscaledpoints = DelaunayMeshes.unscalepoints(mesh, [xc'; yc']')"
},

{
    "location": "UserGuide.html#Constraints-1",
    "page": "User Guide",
    "title": "Constraints",
    "category": "section",
    "text": "Constraints can be added to unconstrained triangulations using addconstraint# create an unconstrained triangulation\nmesh = DelaunayMeshes.Mesh()\nDelaunayMeshes.setboundingbox(mesh, [-15.0, 15.0, -15.0, 15.0])\nseed = rand(UInt32)\nnvertices = 200\npoints = rand(Float64, nvertices, 2)*30. - 15.\npush!(mesh, points)  \n\n# push vertices for two concentric circles\nnringvertices = 40\nθ = linspace(0., 2.*π, nringvertices + 1)[1:end-1]\nrinner = 5.\nrouter = 10.\ninnerRing = rinner * [cos(θ'); sin(θ')]'\nouterRing = router * [cos(θ'); sin(θ')]'\npush!(mesh, innerRing)\npush!(mesh, outerRing)\n\n# add constraints\ninnerConstraintVertexList = [x for x in (nvertices+3+nringvertices):-1:(nvertices+3+1)]\nouterConstraintVertexList = [x for x in (nvertices+3+nringvertices+1):(nvertices + 3 + 2*nringvertices)]\nDelaunayMeshes.addconstraint!(mesh, innerConstraintVertexList)\nDelaunayMeshes.addconstraint!(mesh, outerConstraintVertexList)All faces that are located to the right of the enclosing polygon are marked as being exterior to the triangulated region. The current status of a face is stored in Mesh.faceLocation. We can use it, for example, to plot all Voronoi vertices that are interior to the region:# get Delaunay edges and Voronoi vertices\nxc, yc = DelaunayMeshes.getdelaunaycoordinates(mesh.tesselation)\nvorVert = DelaunayMeshes.getvoronoivertices(mesh)\n\n# plot all edges and all interior Voronoi vertices\nimport Winston\np = Winston.FramedPlot(aspect_ratio=1)\nWinston.add(p, Winston.Curve(xc, yc))\nWinston.add(p, Winston.Points( vorVert[mesh.faceLocation,1], vorVert[mesh.faceLocation,2], kind=\"circle\", color=\"green\"))(Image: Constrained triangulation)"
},

{
    "location": "UserGuide.html#Quad-Edge-data-structure-1",
    "page": "User Guide",
    "title": "Quad-Edge data structure",
    "category": "section",
    "text": "Internally, DelaunayMeshes uses a quad-edge data structure to maintain the triangulation. Edges have triangulation vertices as sources while the source of a dual edge is a face. All basic quad-edge operations are provided and can be used to match the face index to a triangle (and vice versa). The user may consult the API documentation for details."
},

{
    "location": "index.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DelaunayMeshes",
    "page": "API",
    "title": "DelaunayMeshes",
    "category": "Module",
    "text": "This module provides high-level functionality to generate, manage and refine meshes based on constrained Delaunay triangulations.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.addconstraint!-Tuple{DelaunayMeshes.Mesh,Array{Int64,1}}",
    "page": "API",
    "title": "DelaunayMeshes.addconstraint!",
    "category": "Method",
    "text": "addconstraint!(mesh::Mesh, vertices::Vector{Int})\n\nAdds a constraint, given as a list of existing vertices, into the DelaunayMeshes.\n\nRemarks\n\nConstraints are closed polygons. A constraint edge from the last vertex to the first vertex is automatically inserted to ensure closure.\nAll faces interior of the constrained region are marked as \"exterior\" to the DelaunayMeshes. A face is interior to the constrained region if it is located to the right of the given constraint polygon (i.e. the constraint vertices circle the interior region counter-clockwise).\nConstraint edges must not intersect. For performance reasons, this is not checked.\nAny vertex of the triangulation can intersect either zero or two constraint edges.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.getvoronoivertices-Tuple{DelaunayMeshes.Mesh}",
    "page": "API",
    "title": "DelaunayMeshes.getvoronoivertices",
    "category": "Method",
    "text": "getvoronoivertices(mesh::Mesh)\n\nComputes the Voronoi vertices of the given DelaunayMeshes.\n\nRemarks\n\nThe results are cached but need to be re-computed whenever a new vertex or a new constraint is pushed into the DelaunayMeshes.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.ConstraintType",
    "page": "API",
    "title": "DelaunayMeshes.ConstraintType",
    "category": "Type",
    "text": "Implementation of the QuadEdge-Datastructure. [Guibas85]_\n\n.. [Guibas85] Leonidas Guibas and Jorge Stolfi, Primitives for the manipulation of general subdivisions and the computation of Voronoi diagrams, ACM Transactions on Graphics, 4(2), 1985, 75-123.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.DelaunayTesselation-Tuple{}",
    "page": "API",
    "title": "DelaunayMeshes.DelaunayTesselation",
    "category": "Method",
    "text": "DelaunayTesselation()\n\nCreates an initial triangulation consisting of a maximal triangle and the corresponding faces.\n\nRemarks\n\nThe inner face coordinates are the centroid of the triangle.\nThe outer face coordinates are located at ( ).\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.Mesh",
    "page": "API",
    "title": "DelaunayMeshes.Mesh",
    "category": "Type",
    "text": "Defines a mesh. A mesh consists of a bounding box providing bounds to the vertex coordinates, constraints consisting of closed polygons and the corresponding Delaunay tesselation.\n\n\n\n"
},

{
    "location": "index.html#Base.convert-Tuple{Type{DiffEqPDEBase.SimpleFEMMesh},DelaunayMeshes.Mesh}",
    "page": "API",
    "title": "Base.convert",
    "category": "Method",
    "text": "convert(::Type{DiffEqPDEBase.SimpleFEMMesh}, mesh::Mesh)\n\nConverts the given mesh to a format that can be used by the FEM-solvers.\n\n\n\n"
},

{
    "location": "index.html#Base.push!-Tuple{DelaunayMeshes.Mesh,Array{Float64,2}}",
    "page": "API",
    "title": "Base.push!",
    "category": "Method",
    "text": "push!(mesh::Mesh, points::Array{Float64, 2})\n\nInserts a number of vertices into the mesh.\n\nRemarks\n\nIf the mesh was initialized without a bounding box, it is determined by the first set of vertices that are pushed into the triangualation. All consecutive push operations cannot insert vertices that exceed the bounding box limits\nNew faces that are created by a push operation are marked as interior by default, i.e. it is forbidden to push vertices into an exterior region.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.base-Tuple{Int64}",
    "page": "API",
    "title": "DelaunayMeshes.base",
    "category": "Method",
    "text": "base(ei::Int)\n\nComputes base index of index ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.checkconvexityquadriliteral-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Array{Int64,1}}",
    "page": "API",
    "title": "DelaunayMeshes.checkconvexityquadriliteral",
    "category": "Method",
    "text": "checkconvexityquadriliteral(sd::DelaunayTesselation, vertices::Vector{Int})\n\nChecks if the given quadriliteral is strictly convex.\n\nRemarks\n\nThe end points of the quadriliteral must be given in CCW order\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.computescaleandshiftfromboundingbox-Tuple{Array{Float64,1}}",
    "page": "API",
    "title": "DelaunayMeshes.computescaleandshiftfromboundingbox",
    "category": "Method",
    "text": "computescaleandshiftfromboundingbox(boundingBox::Vector{Float64})\n\nComputes the scale and horizontal/vertical shift from a given bounding box. This establishes an aspect-ratio-conserving transformation from the given bounding box into a square that fits comfortably into the outer bounding triangle of the DelaunayMeshes. The square (blue) fills ten per cent of the maximally inscribed square (dashed).\n\n(Image: inscribed)\n\nSpecifically, if (c_x c_y) denote the center of the bounding box and Delta = max((x_mathrmmax - x_mathrmmin) (x_mathrmmax - x_mathrmmin)) is the maximum extent of the bounding box, then the coordinate transformation to the new primed coordinates is given by (x y) = 01 (x - c_x y - c_y)Delta + (c_x c_y), where (c_x c_y) = (15 125) denotes the center of the inscribed square.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.connect!-Tuple{DelaunayMeshes.SubDivision{T},Int64,Int64,T}",
    "page": "API",
    "title": "DelaunayMeshes.connect!",
    "category": "Method",
    "text": "connect!{T}(sd::SubDivision, ai::Int, bi::Int, newFace::T)\n\nConnects the end vertex of edge ai to the start vertex bi.\n\nRemarks\n\nThis operation adds another face to the subdivision. The new face will be to the right of the new edge and will be correctly assigned to all adjacent edges.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.dest-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.dest",
    "category": "Method",
    "text": "dest(sd::SubDivision, ei::Int)\n\nGets the index of the destination vertex of edge index ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.dprev-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.dprev",
    "category": "Method",
    "text": "dprev(sd::SubDivision, ei::Int)\n\nGets the next edge ending at the destination of ei in clockwise direction.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.endpoints!-Tuple{DelaunayMeshes.SubDivision,Int64,Int64,Int64,Int64,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.endpoints!",
    "category": "Method",
    "text": "endpoints!(sd::SubDivision, ei::Int, vio::Int, vid::Int, lfi::Int, rfi::Int)\n\nAssigns the start(end) point vio(vid) and the left(right) faces lfi(rfi) to given edge ei.\n\nRemarks\n\nThe start(end) points of the symmetric and symmetric dual edges are assigned correctly.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.findintersectingedges-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Int64,Int64,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.findintersectingedges",
    "category": "Method",
    "text": "findintersectingedges(sd::DelaunayTesselation, v1::Int, v2::Int, eg::Int)\n\nFinds all edges of the given tesselation that intersect with a virtual edge from v1 to v2.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.getadjacentquadriliteral-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Int64}",
    "page": "API",
    "title": "DelaunayMeshes.getadjacentquadriliteral",
    "category": "Method",
    "text": "getadjacentquadriliteral(sd::DelaunayTesselation, ei::Int)\n\nGets vertices of quadriliteral formed by the two triangles adjacent to given edge in ccw order.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.getalltriangles-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}}}",
    "page": "API",
    "title": "DelaunayMeshes.getalltriangles",
    "category": "Method",
    "text": "getalltriangles(sd::SubDivision)\n\nReturns all triangles as edge lists.\n\nRemarks\n\nThe boundary edges are included.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.getfacesinsideregion-Tuple{DelaunayMeshes.Mesh,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.getfacesinsideregion",
    "category": "Method",
    "text": "getfacesinsideregion(mesh::Mesh, ei::Int)\n\nFinds all faces that are located inside the same region as the face identified by org(ei).\n\nRemarks\n\nThe algorithm performs a breadth-first search over all adjacent triangles while omitting triangles that requires crossing a boundary.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.incircle-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Int64,Int64,Int64,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.incircle",
    "category": "Method",
    "text": "incircle(sd::DelaunayTesselation, v1::Int, v2::Int, v3::Int, v4::Int)\n\nTests if the vertex v4 is strictly inside the circumcircle of the triangle v1-v2-v3.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.insertconstraint!-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Int64,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.insertconstraint!",
    "category": "Method",
    "text": "insertconstraint!(sd::DelaunayTesselation, v1::Int, v2::Int)\n\nInserts constraint going from v1 to v2.\n\nRemarks\n\nThe constrained region is assumed to be on the right of the directed edge.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.insertpoint!-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Array{Float64,1},Int64}",
    "page": "API",
    "title": "DelaunayMeshes.insertpoint!",
    "category": "Method",
    "text": "insertpoint!(sd::DelaunayTesselation, x::Point, ei::Int)\n\nInserts a given point into an existing DelaunayMeshes.\n\nRemarks\n\nThe algorithm identifies the triangle in the existing triangulation that contains the vertex to be inserted, taking <paramref name =\"ei\"/> as start point for the search. It then connects the newly-inserted point to all vertices of the enclosing triangle.\nIf an existing vertex is already located at the point coordinates an error is thrown.\n\n#Returns\n\nnew edge, last edge of target triangle (lnext(lastEdge) = newEdge)\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.invrot-Tuple{Int64}",
    "page": "API",
    "title": "DelaunayMeshes.invrot",
    "category": "Method",
    "text": "invrot(ei::Int)\n\nInverse operation of rot(ei).\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.lface-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.lface",
    "category": "Method",
    "text": "lface(sd::SubDivision, ei::Int)\n\nReturns the index to face left of given edge ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.lnext-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.lnext",
    "category": "Method",
    "text": "lnext(sd::SubDivision, ei::Int)\n\nGets the next edge around the left face of ei and starting at destination of ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.locatevertex-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Array{Float64,1},Int64}",
    "page": "API",
    "title": "DelaunayMeshes.locatevertex",
    "category": "Method",
    "text": "locatevertex(sd::DelaunayTesselation, x::Point, ei::Int)\n\nLocates point x in existing subdivision by traversing the triangulation starting at edge ei.\n\nRemarks\n\nThe point coordinates must be exactly inside the enclosing triangle.\nThe algorithm used here [Brown]_ is more stable than the original version\n\n.. [Brown] Brown, Faigle: A robust efficient algorithm for point location in triangulations\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.lprev-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.lprev",
    "category": "Method",
    "text": "lprev(sd::SubDivision, ei::Int)\n\nGets the previous edge around the left face of ei and ending at source of ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.makeedge!-Tuple{DelaunayMeshes.SubDivision}",
    "page": "API",
    "title": "DelaunayMeshes.makeedge!",
    "category": "Method",
    "text": "makeedge!(sd::SubDivision)\n\nCreates an empty edge and adds it to given SubDivision.\n\n#returns\n\nNumber of newly-added edge.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.onext-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.onext",
    "category": "Method",
    "text": "onext(sd::SubDivision, ei::Int)\n\nGets the next edge starting at the origin of ei in counter-clockwise direction.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.oprev-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.oprev",
    "category": "Method",
    "text": "oprev(sd::SubDivision, ei::Int)\n\nGets the next edge starting at the origin of ei in clockwise direction.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.org-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.org",
    "category": "Method",
    "text": "org(sd::SubDivision, ei::Int)\n\nGets the index of the origin vertex of edge index ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.pointsintersect-Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "API",
    "title": "DelaunayMeshes.pointsintersect",
    "category": "Method",
    "text": "pointsintersect(p1::Point, q1::Point, p2::Point, q2::Point)\n\nChecks if line segments given by their end points intersect.\n\nRemarks\n\nTaken from geeksforgeeks.org/check-if-two-given-line-segments-intersect\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.removedeletededges!-Tuple{DelaunayMeshes.SubDivision}",
    "page": "API",
    "title": "DelaunayMeshes.removedeletededges!",
    "category": "Method",
    "text": "removedeletededges!(sd::SubDivision)\n\nRemoves all edges that were marked as deleted in DelaunayMeshes.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.removeintersectingedges!-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Int64,Int64,Array{Int64,1}}",
    "page": "API",
    "title": "DelaunayMeshes.removeintersectingedges!",
    "category": "Method",
    "text": "removeintersectingedges!(\n  sd::DelaunayTesselation, v1::Int, v2::Int, intersecting::Vector{Int})\n\nRemoves all intersecting edges by iteratively swapping diagonals of enclosing quadriliterals.\n\nRemarks\n\nFor details see Sloan: A fast algorithm for generating constrained Delaunay triangulations.\nReturns the index of the constraint edge\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.restoredelaunay!-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Int64,Int64,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.restoredelaunay!",
    "category": "Method",
    "text": "restoredelaunay!(sd::DelaunayTesselation, ei::Int, si::Int, xi::Int)\n\nRestores Delaunay constraint for edge structure by swapping edges of adjacent triangles.\n\nRemarks\n\nIt examines all triangles to check if the new vertex xi is inside its circumcircle.\nThe algorithm cycles through all adjacent edges starting at startEdge and stopping at stopEdge.\nIf an edge is swapped, its neighbours are examined as well.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.rface-Tuple{DelaunayMeshes.SubDivision,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.rface",
    "category": "Method",
    "text": "lface(sd::SubDivision, ei::Int)\n\nReturns the index to face left of given edge ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.rightof-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Array{Float64,1},Int64}",
    "page": "API",
    "title": "DelaunayMeshes.rightof",
    "category": "Method",
    "text": "rightof(sd::DelaunayTesselation, vi::Int, ei::Int)\nrightof(sd::DelaunayTesselation, x::Point, ei::Int)\n\nComputes the location of vertex vi relative to line defined by ei.\n\nRemarks\n\nInternally, the computation is forwarded to GeometricalPredicates.orientation.\n\nReturns\n\nTrue, if vi is located to the right of the edge\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.rot-Tuple{Int64}",
    "page": "API",
    "title": "DelaunayMeshes.rot",
    "category": "Method",
    "text": "rot(ei::Int)\n\nComputes the index of the dual edge to index ei.\n\nRemarks\n\nThe dual edge is obtained by rotating ei counter-clockwise by 90 degrees.\nThe start(end) points of the dual edge can be interpreted as representing the right(left) face.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.splice!-Tuple{DelaunayMeshes.SubDivision,Int64,Int64}",
    "page": "API",
    "title": "DelaunayMeshes.splice!",
    "category": "Method",
    "text": "splice!(sd::SubDivision, ai::Int, bi::Int)\n\nCombines or splices two given edge rings ai and bi.\n\nRemarks\n\nIf the origins of ai and bi are identical, it cuts the two edge rings\nIf the origins differ, it joins the two edge rings\nAfter the operation, the new edge rings satisfy the identities onext(ai) = onext(bi) and onext(bi) = onext(ai) and the corresponding identities for the dual rings.\nThe origins of the edges are neither evaluated nor are theyaffected by this operation and need to be updated manually.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.swap!-Tuple{DelaunayMeshes.SubDivision{Array{Float64,1}},Int64}",
    "page": "API",
    "title": "DelaunayMeshes.swap!",
    "category": "Method",
    "text": "swap!(sd::DelaunayTesselation, ei::Int)\n\nSwaps edge ei such that the new edge connects the apexes of the triangles adjacent to the old edge.\n\nRemarks\n\nSpecifially, if a = oprev(ei) and b = oprev(sym(ei)) before the swap, then the edge ei will connect dest(a) with dest(b) after the swap.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes.sym-Tuple{Int64}",
    "page": "API",
    "title": "DelaunayMeshes.sym",
    "category": "Method",
    "text": "sym(ei::Int)\n\nComputes index of symmetric index ei.\n\n\n\n"
},

{
    "location": "index.html#DelaunayMeshes-1",
    "page": "API",
    "title": "DelaunayMeshes",
    "category": "section",
    "text": "Modules = [DelaunayMeshes]"
},

]}
