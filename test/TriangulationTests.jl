import GeometricalPredicates
#import Winston

# turn off matplotlib's gui
#PyPlot.pygui(false)

@testset "TriangulationTests" begin
  @testset "Iterator" begin
    local tess =
     setuptwopointtesselation(
       [1.375, 1.25],
       [1.5,1.5])

    local faces = Array{DelaunayMeshes.VertexIndex, 1}()
    for face in tess
      push!(faces, face.face)
    end

    @test sort(faces) == [1,2,3,4,5,6]
  end

  @testset "InsertVertexOnEdge" begin
    local tess =
     setuptwopointtesselation(
       [1.375, 1.25],
       [1.5,1.5])

    # push vertex on edge from (1.,1.) to (1.375, 1.25)
    local e14 = findedgeconnectingvertices(tess, 1, 4)
    @test DelaunayMeshes.dest(tess, e14) == 4

    local ev = 0.5*(tess.vertices[1] + tess.vertices[4])
    push!(tess, ev')

    DelaunayMeshes.removedeletededges!(tess)

    # check that there is no edge from 1->4 but there are connecting edges
    @test_throws AssertionError findedgeconnectingvertices(tess, 1, 4)

    local e16 = findedgeconnectingvertices(tess, 1, 6)
    @test DelaunayMeshes.dest(tess, e16) == 6

    local e64 = findedgeconnectingvertices(tess, 6, 4)
    @test DelaunayMeshes.dest(tess, e64) == 4

    local e62 = findedgeconnectingvertices(tess, 6, 2)
    @test DelaunayMeshes.dest(tess, e62) == 2

    local e65 = findedgeconnectingvertices(tess, 6, 5)
    @test DelaunayMeshes.dest(tess, e65) == 5

    #=
    # plot it
    xc, yc = getplotedges(tess)
    local vertexPoints = reduce(hcat, tess.vertices)'
    p = Winston.FramedPlot(aspect_ratio=1
      #,xrange=[1.496 1.50], yrange=[1.256, 1.26]
      )
    Winston.add(p, Winston.Points(vertexPoints[:,1], vertexPoints[:,2], kind="circle", color="red"))
    Winston.add(p, Winston.Curve(xc, yc))
    Winston.savefig(p, "InsertVertexOnEdge.svg")
    =#
  end

  @testset "TriangulateWithConstraints" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    # seed rng and store seed in file
    #local seed = 3393535361
    local seed = rand(UInt32)
    local fs = open("TriangulateWithConstraints.seed", "w")
    write(fs, string(seed))
    close(fs)

    # push some random vertices
    srand(seed)
    local nvertices = 20
    local points = rand(Float64, nvertices, 2)*0.1 + 1.45
    push!(tess, points)

    # push constraint vertices
    local constraintVertices = [
      -1.0 -1.0;
      -0.5 -1.0;
      -0.5 0.0;
      0.5 0.0;
      0.5 -1.0;
      1.0 -1.0;
      1.0 1.0;
      -1.0 1.0;
    ]
    local scaledConstrVert = 1.5*ones(constraintVertices) + 0.025*constraintVertices
    push!(tess, scaledConstrVert)
    @test length(tess.vertices) == (nvertices + 8 + 3)

    local faulty = testdelaunayness(tess)
    @test length(faulty) == 0

    # constraints
    local n = nvertices + 3

    local constraints = [
      n+1 n+2;
      n+2 n+3;
      n+3 n+4;
      n+4 n+5;
      n+5 n+6;
      n+6 n+7;
      n+7 n+8;
      n+8 n+1;
    ]
    for i=1:size(constraints, 1)
      DelaunayMeshes.DelaunayMeshes.insertconstraint!(tess, constraints[i, 1], constraints[i,2])
    end

    #=
    # plot it
    xc, yc = DelaunayMeshes.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("TriangulateWithConstraints.svg")
    =#

    validateconstraints(tess, constraints)
    testvertexcache(tess)
  end

  @testset "DelaunayMeshes.swap" begin
    local tess =
     setuptwopointtesselation(
       convert(DelaunayMeshes.Point, [1.375, 1.25]),
       convert(DelaunayMeshes.Point, [1.5,1.5]))

    local ea = 1
    local eb = 5
    local ec = 9

    # test initial triangulation connectivity
    for e in  DelaunayMeshes.lface_edges(tess, ea)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 1
    end
    for e in DelaunayMeshes.lface_edges(tess, 19)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 3
    end
    for e in DelaunayMeshes.lface_edges(tess, 13)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 4
    end
    for e in DelaunayMeshes.lface_edges(tess, 5)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 5
    end
    for e in DelaunayMeshes.lface_edges(tess, 25)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 6
    end

    # swap edge
    local ei = 21
    DelaunayMeshes.swap!(tess, ei)

    # test new connectivity
    for e in DelaunayMeshes.lface_edges(tess, ea)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 1
    end
    for e in DelaunayMeshes.lface_edges(tess, 19)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 3
    end
    for e in DelaunayMeshes.lface_edges(tess, DelaunayMeshes.sym(ei))
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 4
    end
    for e in DelaunayMeshes.lface_edges(tess, 5)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 5
    end
    for e in DelaunayMeshes.lface_edges(tess, ei)
      local lf = DelaunayMeshes.lface(tess, e)
      @test lf == 6
    end
  end

  @testset "rightof" begin
    local tess =
     setuptwopointtesselation(
       [1.375, 1.25],
       [1.5,1.5])

    local ea = 1
    local eb = 5
    local ec = 9

    @test DelaunayMeshes.rightof(tess, 4, ea) == false
    @test DelaunayMeshes.rightof(tess, 4, DelaunayMeshes.sym(ea)) == true
    @test DelaunayMeshes.rightof(tess, 4, eb) == false
    @test DelaunayMeshes.rightof(tess, 4, DelaunayMeshes.sym(eb)) == true
    @test DelaunayMeshes.rightof(tess, 4, ec) == false
    @test DelaunayMeshes.rightof(tess, 4, DelaunayMeshes.sym(ec)) == true
  end

  @testset "locatevertex" begin
    local tess =
     setuptwopointtesselation(
       [1.375, 1.25],
       [1.5,1.5])

    local x = [1.5, 1.25]
    local y = [1.65, 1.5]

    @test locateandassertintriangle(tess, x, 1)
    @test locateandassertintriangle(tess, x, 33)
    @test locateandassertintriangle(tess, x, 21)
    @test locateandassertintriangle(tess, y, 1)
    @test locateandassertintriangle(tess, y, 33)
    @test locateandassertintriangle(tess, y, 21)

    # try to locate existing vertex
    local e1 = DelaunayMeshes.locatevertex(tess, [1.375, 1.25], 1)
    @test (DelaunayMeshes.org(tess, e1) == 4 || DelaunayMeshes.dest(tess, e1) == 4)
    local e2 = DelaunayMeshes.locatevertex(tess, [1.375, 1.25], 5)
    @test (DelaunayMeshes.org(tess, e2) == 4 || DelaunayMeshes.dest(tess, e2) == 4)
    local e3 = DelaunayMeshes.locatevertex(tess, [1.375, 1.25], 7)
    @test (DelaunayMeshes.org(tess, e3) == 4 || DelaunayMeshes.dest(tess, e3) == 4)
    local e4 = DelaunayMeshes.locatevertex(tess, [1.375, 1.25], 33)
    @test (DelaunayMeshes.org(tess, e4) == 4 || DelaunayMeshes.dest(tess, e4) == 4)
    local e5 = DelaunayMeshes.locatevertex(tess, [1.375, 1.25], 35)
    @test (DelaunayMeshes.org(tess, e5) == 4 || DelaunayMeshes.dest(tess, e5) == 4)
  end

  @testset "onlinesegment" begin
    local tess =
     setuptwopointtesselation(
       [1.375, 1.25],
       [1.5,1.5])

    push!(tess.vertices, [1.75, 1.75])
    push!(tess.vertices, [1.25, 1.25])
    push!(tess.vertices, [1.0, 1.0] + 2.0*[0.375, 0.25])
    push!(tess.vertices, [1.0, 1.0] + 0.75*[0.375, 0.25])

    @test DelaunayMeshes.onlinesegment(tess, 6, 1, 5) == false
    @test DelaunayMeshes.onlinesegment(tess, 7, 1, 5) == true
    @test DelaunayMeshes.onlinesegment(tess, 1, 1, 5) == true
    @test DelaunayMeshes.onlinesegment(tess, 5, 1, 5) == true

    @test DelaunayMeshes.onlinesegment(tess, 8, 1, 4) == false
    @test DelaunayMeshes.onlinesegment(tess, 9, 1, 4) == true
    @test DelaunayMeshes.onlinesegment(tess, 1, 1, 4) == true
    @test DelaunayMeshes.onlinesegment(tess, 4, 1, 4) == true
  end

  @testset "DeleteEdge" begin
    local tess =
     setuptwopointtesselation(
       [1.375, 1.25],
       [1.5,1.5])

    #=
    # plot it
    xc, yc = getplotedges(tess)
    local vertexPoints = reduce(hcat, tess.vertices)'
    p = Winston.FramedPlot(aspect_ratio=1
      #,xrange=[1.496 1.50], yrange=[1.256, 1.26]
      )
    Winston.add(p, Winston.Points(vertexPoints[:,1], vertexPoints[:,2], kind="circle", color="red"))
    Winston.add(p, Winston.Curve(xc, yc))
    Winston.savefig(p, "DeleteEdge.svg")
    =#

    local e25 = findedgeconnectingvertices(tess, 4, 5)
    DelaunayMeshes.deleteedge!(tess, e25)
    DelaunayMeshes.removedeletededges!(tess)
    testvertexcache(tess)

    local e19 = findedgeconnectingvertices(tess, 4, 2)
    local e23 = findedgeconnectingvertices(tess, 4, 3)
    local e29 = findedgeconnectingvertices(tess, 2, 5)
    local n19 = DelaunayMeshes.onext(tess, e19)
    local p23 = DelaunayMeshes.oprev(tess, e23)
    local n29 = DelaunayMeshes.lnext(tess, e29)

    @test DelaunayMeshes.org(tess, n19) == 4
    @test DelaunayMeshes.dest(tess, n19) == 3
    @test DelaunayMeshes.org(tess, p23) == 4
    @test DelaunayMeshes.dest(tess, p23) == 2
    @test DelaunayMeshes.org(tess, n29) == 5
    @test DelaunayMeshes.dest(tess, n29) == 3
    @test n19 == e23
    @test p23 == e19
    @test_throws AssertionError findedgeconnectingvertices(tess, 4, 5)

    # check that face is correctly assigned
    local faceEdges = DelaunayMeshes.lface_edges(tess, 21)
    for curEdge in faceEdges
      @test DelaunayMeshes.lface(tess, curEdge) == 6
    end
  end

  @testset "InsertVertex" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    # insert a vertex and check that it is connected to the boundary vertices
    local x = [1.5, 1.5]
    DelaunayMeshes.insertpoint!(tess, x, 1)

    @test tess.vertices[4] ≈ x
    local e1 = tess.vertexCache[4]
    local ed = DelaunayMeshes.dest(tess, e1)
    @test DelaunayMeshes.org(tess, e1) == 4
    @test ed ∈ [1,2,3]
    local e1lf = DelaunayMeshes.lface(tess, e1)
    local e1rf = DelaunayMeshes.rface(tess, e1)
    @test e1lf != 2
    @test e1rf != 2
    @test e1lf != e1rf

    local e2 = DelaunayMeshes.onext(tess, e1)
    ed = ed == 3 ? 1 : ed + 1
    @test e2 != e1
    @test DelaunayMeshes.org(tess, e2) == 4
    @test DelaunayMeshes.dest(tess, e2) == ed
    local e2lf = DelaunayMeshes.lface(tess, e2)
    local e2rf = DelaunayMeshes.rface(tess, e2)
    @test e1lf == e2rf
    @test e2lf != e1lf

    local e3 = DelaunayMeshes.onext(tess, e2)
    ed = ed == 3 ? 1 : ed + 1
    @test e3 != e1
    @test e3 != e2
    @test DelaunayMeshes.org(tess, e3) == 4
    @test DelaunayMeshes.dest(tess, e3) == ed
    local e3lf = DelaunayMeshes.lface(tess, e3)
    local e3rf = DelaunayMeshes.rface(tess, e3)
    @test e2lf == e3rf
    @test e3lf == e1rf

    @test DelaunayMeshes.onext(tess, e3) == e1
  end

  @testset "InsertExistingVertex" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    # insert a vertex
    local x = [1.5, 1.5]
    DelaunayMeshes.insertpoint!(tess, x, 1)

    # try to insert same vertex again
    @test_throws ErrorException DelaunayMeshes.insertpoint!(tess, x, 1)
  end

  @testset "DelaunayMeshes.restoredelaunay" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    # insert a triangle and a point inside its circumcircle
    local x1 = [1.5, 1.5] + 0.2 * [cos(0), sin(0)]
    local x2 = [1.5, 1.5] + 0.2 * [cos(pi), sin(pi)]
    local x3 = [1.5, 1.5] + 0.2 * [cos(pi/2), sin(pi/2)]
    local x4 = [1.5, 1.5] + 0.19 * [cos(pi/4), sin(pi/4)]

    DelaunayMeshes.insertpoint!(tess, x1, 1)
    DelaunayMeshes.insertpoint!(tess, x2, 1)
    DelaunayMeshes.insertpoint!(tess, x3, 1)
    stopEi, startEi = DelaunayMeshes.insertpoint!(tess, x4, 1)

    DelaunayMeshes.DelaunayMeshes.restoredelaunay!(tess, startEi, stopEi, 7)

    local ostar = DelaunayMeshes.ostar(tess, 5)
    local ostarvert = map((ei) -> DelaunayMeshes.dest(tess, ei), ostar)
    @test sort(ostarvert) == [1,3,4,6,7]

    # plot it
    #=
    xc, yc = DelaunayMeshes.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("before_swap.svg")
    =#

  end

  @testset "Push" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    # seed rng and store seed in file
    local seed = rand(UInt32)
    local fs = open("push.seed", "w")
    write(fs, string(seed))
    close(fs)

    local nvertices = 20
    local points = rand(Float64, nvertices, 2)*1e-4 + 1.5
    push!(tess, points)

    #=
    # plot it
    xc, yc = DelaunayMeshes.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("random_triangles.svg")
    =#

    local faulty = testdelaunayness(tess)
    @test length(faulty) == 0
  end

  @testset "Intersects" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    # fill in a couple of vertices
    push!(tess.vertices, [1.4, 1.4])
    push!(tess.vertices, [1.6, 1.6])
    push!(tess.vertices, [1.4, 1.6])
    push!(tess.vertices, [1.6, 1.4])
    push!(tess.vertices, [1.525, 1.525])
    push!(tess.vertices, [1.65, 1.65])
    push!(tess.vertices, [1.5, 1.5])
    push!(tess.vertices, [1.5, 1.45])
    push!(tess.vertices, [1.5, 1.425])
    push!(tess.vertices, [1.45, 1.49])

    push!(tess.vertexCache, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)

    # and insert some test edges
    local e45 = insertsingleedge!(tess, 4, 5, 1, 1);
    local e67 = insertsingleedge!(tess, 6, 7, 1, 1);
    local e89 = insertsingleedge!(tess, 8, 9, 1, 1);
    local e1011 = insertsingleedge!(tess, 10, 11, 1, 1);
    local e1213 = insertsingleedge!(tess, 12, 13, 1, 1);

    #=
    # plot it
    xc, yc = DelaunayMeshes.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("intersecting.svg")
    =#

    # intersect
    @test DelaunayMeshes.edgeintersects(tess, e45, e67)
    @test DelaunayMeshes.edgeintersects(tess, e67, e45)
    @test DelaunayMeshes.edgeintersects(tess, DelaunayMeshes.sym(e45), e67)
    @test DelaunayMeshes.edgeintersects(tess, e67, DelaunayMeshes.sym(e45))
    @test DelaunayMeshes.edgeintersects(tess, e45, DelaunayMeshes.sym(e67))
    @test DelaunayMeshes.edgeintersects(tess, DelaunayMeshes.sym(e67), e45)
    @test DelaunayMeshes.edgeintersects(tess, DelaunayMeshes.sym(e67), DelaunayMeshes.sym(e45))
    @test DelaunayMeshes.edgeintersects(tess, DelaunayMeshes.sym(e45), DelaunayMeshes.sym(e67))

    # coplanar intersect
    @test DelaunayMeshes.edgeintersects(tess, e89, e45)
    @test DelaunayMeshes.edgeintersects(tess, e45, e89)

    # endpoints intersect
    @test DelaunayMeshes.edgeintersects(tess, e1011, e45)

    # don't intersect
    @test !DelaunayMeshes.edgeintersects(tess, e1011, e1213)
    @test !DelaunayMeshes.edgeintersects(tess, e1213, e67)

  end

  @testset "findintersectingedges" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    local points = [1.5 1.5; 1.51 1.45; 1.49 1.425; 1.45 1.49; 1.525 1.5; 1.48 1.48; 1.42 1.44]
    push!(tess, points)
    testvertexcache(tess)

    #=
    # plot it
    xc, yc = DelaunayMeshes.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("findintersecting.svg")
    =#

    # check some intersections for non-existing edges
    intersecting, alreadyIn = DelaunayMeshes.findintersectingedges(tess, 6, 8, 1)
    @test !alreadyIn
    @test length(intersecting) == 2
    @test edgeconnectsvertices(tess, intersecting[1], 5, 9)
    @test edgeconnectsvertices(tess, intersecting[2], 5, 4)

    intersecting, alreadyIn = DelaunayMeshes.findintersectingedges(tess, 10, 8, 1)
    @test !alreadyIn
    @test length(intersecting) == 3
    @test edgeconnectsvertices(tess, intersecting[1], 6, 9)
    @test edgeconnectsvertices(tess, intersecting[2], 5, 9)
    @test edgeconnectsvertices(tess, intersecting[3], 5, 4)

    # check some edges that are already in
    intersecting, alreadyIn = DelaunayMeshes.findintersectingedges(tess, 6, 9, 1)
    @test alreadyIn
    @test length(intersecting) == 1
    @test edgeconnectsvertices(tess, intersecting[1], 6, 9)

    intersecting, alreadyIn = DelaunayMeshes.findintersectingedges(tess, 9, 5, 1)
    @test alreadyIn
    @test length(intersecting) == 1
    @test edgeconnectsvertices(tess, intersecting[1], 9, 5)
  end

  @testset "RemoveIntersectingEdges" begin
    local tess = DelaunayMeshes.DelaunayTesselation()

    local points = [1.5 1.5; 1.51 1.45; 1.49 1.425; 1.45 1.49; 1.525 1.5; 1.48 1.48; 1.42 1.44]
    push!(tess, points)
    testvertexcache(tess)

    intersecting, alreadyIn = DelaunayMeshes.findintersectingedges(tess, 10, 8, 1)
    @test !alreadyIn
    @test length(intersecting) == 3

    local constraint = DelaunayMeshes.removeintersectingedges!(tess, 10, 8, intersecting)
    testvertexcache(tess)
    @test length(intersecting) == 3 # it should not affect the container

    # check that constraint connects the correct vertices and does not
    # intersect any other edges
    @test edgeconnectsvertices(tess, constraint, 10, 8)
    for i in 1:4:length(tess.edges)
      if (DelaunayMeshes.base(i) != DelaunayMeshes.base(constraint)) && (DelaunayMeshes.org(tess, i) != DelaunayMeshes.org(tess, constraint)) && (DelaunayMeshes.org(tess, i) != DelaunayMeshes.dest(tess, constraint)) && (DelaunayMeshes.dest(tess, i) != DelaunayMeshes.org(tess, constraint)) && (DelaunayMeshes.dest(tess, i) != DelaunayMeshes.dest(tess, constraint))
          @test !DelaunayMeshes.edgeintersects(tess, constraint, i)
      end
    end

    #=
    # plot it
    xc, yc = DelaunayMeshes.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("removeintersecting.svg")
    =#
  end
end
