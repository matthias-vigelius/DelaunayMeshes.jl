using Base.Test

import Triangulation
import GeometricalPredicates
import TestHelpers

import Winston

# turn off matplotlib's gui
#PyPlot.pygui(false)

@testset "TriangulationTests" begin
  @testset "InsertVertexOnEdge" begin
    local tess =
     TestHelpers.SetUpTwoPointTesselation(
       [1.375, 1.25],
       [1.5,1.5])

    # push vertex on edge from (1.,1.) to (1.375, 1.25)
    local e14 = TestHelpers.FindEdgeConnectingVertices(tess, 1, 4)
    @test QuadEdge.Dest(tess, e14) == 4

    local ev = 0.5*(tess.vertices[1] + tess.vertices[4])
    push!(tess, ev')

    QuadEdge.RemoveDeletedEdges!(tess)

    # check that there is no edge from 1->4 but there are connecting edges
    @test_throws AssertionError TestHelpers.FindEdgeConnectingVertices(tess, 1, 4)

    local e16 = TestHelpers.FindEdgeConnectingVertices(tess, 1, 6)
    @test QuadEdge.Dest(tess, e16) == 6

    local e64 = TestHelpers.FindEdgeConnectingVertices(tess, 6, 4)
    @test QuadEdge.Dest(tess, e64) == 4

    local e62 = TestHelpers.FindEdgeConnectingVertices(tess, 6, 2)
    @test QuadEdge.Dest(tess, e62) == 2

    local e65 = TestHelpers.FindEdgeConnectingVertices(tess, 6, 5)
    @test QuadEdge.Dest(tess, e65) == 5

    #=
    # plot it
    xc, yc = TestHelpers.GetPlotEdges(tess)
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
    local tess = Triangulation.DelaunayTesselation()

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

    local faulty = TestHelpers.TestDelaunayness(tess)
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
      Triangulation.InsertConstraint!(tess, constraints[i, 1], constraints[i,2])
    end

    #=
    # plot it
    xc, yc = Triangulation.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("TriangulateWithConstraints.svg")
    =#

    TestHelpers.ValidateConstraints(tess, constraints)
    TestHelpers.TestVertexCache(tess)
  end

  @testset "Swap" begin
    local tess =
     TestHelpers.SetUpTwoPointTesselation(
       convert(Triangulation.Point, [1.375, 1.25]),
       convert(Triangulation.Point, [1.5,1.5]))

    local ea = 1
    local eb = 5
    local ec = 9

    # test initial triangulation connectivity
    for e in  QuadEdge.LFaceEdges(tess, ea)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 1
    end
    for e in QuadEdge.LFaceEdges(tess, 19)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 3
    end
    for e in QuadEdge.LFaceEdges(tess, 13)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 4
    end
    for e in QuadEdge.LFaceEdges(tess, 5)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 5
    end
    for e in QuadEdge.LFaceEdges(tess, 25)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 6
    end

    # swap edge
    local ei = 21
    Triangulation.Swap!(tess, ei)

    # test new connectivity
    for e in QuadEdge.LFaceEdges(tess, ea)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 1
    end
    for e in QuadEdge.LFaceEdges(tess, 19)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 3
    end
    for e in QuadEdge.LFaceEdges(tess, QuadEdge.Sym(ei))
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 4
    end
    for e in QuadEdge.LFaceEdges(tess, 5)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 5
    end
    for e in QuadEdge.LFaceEdges(tess, ei)
      local lf = QuadEdge.LFace(tess, e)
      @test lf == 6
    end
  end

  @testset "RightOf" begin
    local tess =
     TestHelpers.SetUpTwoPointTesselation(
       [1.375, 1.25],
       [1.5,1.5])

    local ea = 1
    local eb = 5
    local ec = 9

    @test Triangulation.RightOf(tess, 4, ea) == false
    @test Triangulation.RightOf(tess, 4, QuadEdge.Sym(ea)) == true
    @test Triangulation.RightOf(tess, 4, eb) == false
    @test Triangulation.RightOf(tess, 4, QuadEdge.Sym(eb)) == true
    @test Triangulation.RightOf(tess, 4, ec) == false
    @test Triangulation.RightOf(tess, 4, QuadEdge.Sym(ec)) == true
  end

  @testset "LocateVertex" begin
    local tess =
     TestHelpers.SetUpTwoPointTesselation(
       [1.375, 1.25],
       [1.5,1.5])

    local x = [1.5, 1.25]
    local y = [1.65, 1.5]

    @test TestHelpers.LocateAndAssertInTriangle(tess, x, 1)
    @test TestHelpers.LocateAndAssertInTriangle(tess, x, 33)
    @test TestHelpers.LocateAndAssertInTriangle(tess, x, 21)
    @test TestHelpers.LocateAndAssertInTriangle(tess, y, 1)
    @test TestHelpers.LocateAndAssertInTriangle(tess, y, 33)
    @test TestHelpers.LocateAndAssertInTriangle(tess, y, 21)

    # try to locate existing vertex
    local e1 = Triangulation.LocateVertex(tess, [1.375, 1.25], 1)
    @test (QuadEdge.Org(tess, e1) == 4 || QuadEdge.Dest(tess, e1) == 4)
    local e2 = Triangulation.LocateVertex(tess, [1.375, 1.25], 5)
    @test (QuadEdge.Org(tess, e2) == 4 || QuadEdge.Dest(tess, e2) == 4)
    local e3 = Triangulation.LocateVertex(tess, [1.375, 1.25], 7)
    @test (QuadEdge.Org(tess, e3) == 4 || QuadEdge.Dest(tess, e3) == 4)
    local e4 = Triangulation.LocateVertex(tess, [1.375, 1.25], 33)
    @test (QuadEdge.Org(tess, e4) == 4 || QuadEdge.Dest(tess, e4) == 4)
    local e5 = Triangulation.LocateVertex(tess, [1.375, 1.25], 35)
    @test (QuadEdge.Org(tess, e5) == 4 || QuadEdge.Dest(tess, e5) == 4)
  end

  @testset "OnLineSegment" begin
    local tess =
     TestHelpers.SetUpTwoPointTesselation(
       [1.375, 1.25],
       [1.5,1.5])

    push!(tess.vertices, [1.75, 1.75])
    push!(tess.vertices, [1.25, 1.25])
    push!(tess.vertices, [1.0, 1.0] + 2.0*[0.375, 0.25])
    push!(tess.vertices, [1.0, 1.0] + 0.75*[0.375, 0.25])

    @test Triangulation.OnLineSegment(tess, 6, 1, 5) == false
    @test Triangulation.OnLineSegment(tess, 7, 1, 5) == true
    @test Triangulation.OnLineSegment(tess, 1, 1, 5) == true
    @test Triangulation.OnLineSegment(tess, 5, 1, 5) == true

    @test Triangulation.OnLineSegment(tess, 8, 1, 4) == false
    @test Triangulation.OnLineSegment(tess, 9, 1, 4) == true
    @test Triangulation.OnLineSegment(tess, 1, 1, 4) == true
    @test Triangulation.OnLineSegment(tess, 4, 1, 4) == true
  end

  @testset "DeleteEdge" begin
    local tess =
     TestHelpers.SetUpTwoPointTesselation(
       [1.375, 1.25],
       [1.5,1.5])

    # plot it
    xc, yc = TestHelpers.GetPlotEdges(tess)
    local vertexPoints = reduce(hcat, tess.vertices)'
    p = Winston.FramedPlot(aspect_ratio=1
      #,xrange=[1.496 1.50], yrange=[1.256, 1.26]
      )
    Winston.add(p, Winston.Points(vertexPoints[:,1], vertexPoints[:,2], kind="circle", color="red"))
    Winston.add(p, Winston.Curve(xc, yc))
    Winston.savefig(p, "DeleteEdge.svg")

    local e25 = TestHelpers.FindEdgeConnectingVertices(tess, 4, 5)
    QuadEdge.DeleteEdge!(tess, e25)
    QuadEdge.RemoveDeletedEdges!(tess)
    TestHelpers.TestVertexCache(tess)

    local e19 = TestHelpers.FindEdgeConnectingVertices(tess, 4, 2)
    local e23 = TestHelpers.FindEdgeConnectingVertices(tess, 4, 3)
    local e29 = TestHelpers.FindEdgeConnectingVertices(tess, 2, 5)
    local n19 = QuadEdge.ONext(tess, e19)
    local p23 = QuadEdge.OPrev(tess, e23)
    local n29 = QuadEdge.LNext(tess, e29)

    @test QuadEdge.Org(tess, n19) == 4
    @test QuadEdge.Dest(tess, n19) == 3
    @test QuadEdge.Org(tess, p23) == 4
    @test QuadEdge.Dest(tess, p23) == 2
    @test QuadEdge.Org(tess, n29) == 5
    @test QuadEdge.Dest(tess, n29) == 3
    @test n19 == e23
    @test p23 == e19
    @test_throws AssertionError TestHelpers.FindEdgeConnectingVertices(tess, 4, 5)

    # check that face is correctly assigned
    local faceEdges = QuadEdge.LFaceEdges(tess, 21)
    for curEdge in faceEdges
      @test QuadEdge.LFace(tess, curEdge) == 6
    end
  end

  @testset "InsertVertex" begin
    local tess = Triangulation.DelaunayTesselation()

    # insert a vertex and check that it is connected to the boundary vertices
    local x = [1.5, 1.5]
    Triangulation.InsertPoint!(tess, x, 1)

    @test tess.vertices[4] ≈ x
    local e1 = tess.vertexCache[4]
    local ed = QuadEdge.Dest(tess, e1)
    @test QuadEdge.Org(tess, e1) == 4
    @test ed ∈ [1,2,3]
    local e1lf = QuadEdge.LFace(tess, e1)
    local e1rf = QuadEdge.RFace(tess, e1)
    @test e1lf != 2
    @test e1rf != 2
    @test e1lf != e1rf

    local e2 = QuadEdge.ONext(tess, e1)
    ed = ed == 3 ? 1 : ed + 1
    @test e2 != e1
    @test QuadEdge.Org(tess, e2) == 4
    @test QuadEdge.Dest(tess, e2) == ed
    local e2lf = QuadEdge.LFace(tess, e2)
    local e2rf = QuadEdge.RFace(tess, e2)
    @test e1lf == e2rf
    @test e2lf != e1lf

    local e3 = QuadEdge.ONext(tess, e2)
    ed = ed == 3 ? 1 : ed + 1
    @test e3 != e1
    @test e3 != e2
    @test QuadEdge.Org(tess, e3) == 4
    @test QuadEdge.Dest(tess, e3) == ed
    local e3lf = QuadEdge.LFace(tess, e3)
    local e3rf = QuadEdge.RFace(tess, e3)
    @test e2lf == e3rf
    @test e3lf == e1rf

    @test QuadEdge.ONext(tess, e3) == e1
  end

  @testset "InsertExistingVertex" begin
    local tess = Triangulation.DelaunayTesselation()

    # insert a vertex
    local x = [1.5, 1.5]
    Triangulation.InsertPoint!(tess, x, 1)

    # try to insert same vertex again
    @test_throws ErrorException Triangulation.InsertPoint!(tess, x, 1)
  end

  @testset "RestoreDelaunay" begin
    local tess = Triangulation.DelaunayTesselation()

    # insert a triangle and a point inside its circumcircle
    local x1 = [1.5, 1.5] + 0.2 * [cos(0), sin(0)]
    local x2 = [1.5, 1.5] + 0.2 * [cos(pi), sin(pi)]
    local x3 = [1.5, 1.5] + 0.2 * [cos(pi/2), sin(pi/2)]
    local x4 = [1.5, 1.5] + 0.19 * [cos(pi/4), sin(pi/4)]

    Triangulation.InsertPoint!(tess, x1, 1)
    Triangulation.InsertPoint!(tess, x2, 1)
    Triangulation.InsertPoint!(tess, x3, 1)
    stopEi, startEi = Triangulation.InsertPoint!(tess, x4, 1)

    Triangulation.RestoreDelaunay!(tess, startEi, stopEi, 7)

    local ostar = QuadEdge.OStar(tess, 5)
    local ostarvert = map((ei) -> QuadEdge.Dest(tess, ei), ostar)
    @test sort(ostarvert) == [1,3,4,6,7]

    # plot it
    #=
    xc, yc = Triangulation.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("before_swap.svg")
    =#

  end

  @testset "Push" begin
    local tess = Triangulation.DelaunayTesselation()

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
    xc, yc = Triangulation.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("random_triangles.svg")
    =#

    local faulty = TestHelpers.TestDelaunayness(tess)
    @test length(faulty) == 0
  end

  @testset "Intersects" begin
    local tess = Triangulation.DelaunayTesselation()

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
    local e45 = TestHelpers.InsertSingleEdge!(tess, 4, 5, 1, 1);
    local e67 = TestHelpers.InsertSingleEdge!(tess, 6, 7, 1, 1);
    local e89 = TestHelpers.InsertSingleEdge!(tess, 8, 9, 1, 1);
    local e1011 = TestHelpers.InsertSingleEdge!(tess, 10, 11, 1, 1);
    local e1213 = TestHelpers.InsertSingleEdge!(tess, 12, 13, 1, 1);

    #=
    # plot it
    xc, yc = Triangulation.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("intersecting.svg")
    =#

    # intersect
    @test Triangulation.EdgeIntersects(tess, e45, e67)
    @test Triangulation.EdgeIntersects(tess, e67, e45)
    @test Triangulation.EdgeIntersects(tess, QuadEdge.Sym(e45), e67)
    @test Triangulation.EdgeIntersects(tess, e67, QuadEdge.Sym(e45))
    @test Triangulation.EdgeIntersects(tess, e45, QuadEdge.Sym(e67))
    @test Triangulation.EdgeIntersects(tess, QuadEdge.Sym(e67), e45)
    @test Triangulation.EdgeIntersects(tess, QuadEdge.Sym(e67), QuadEdge.Sym(e45))
    @test Triangulation.EdgeIntersects(tess, QuadEdge.Sym(e45), QuadEdge.Sym(e67))

    # coplanar intersect
    @test Triangulation.EdgeIntersects(tess, e89, e45)
    @test Triangulation.EdgeIntersects(tess, e45, e89)

    # endpoints intersect
    @test Triangulation.EdgeIntersects(tess, e1011, e45)

    # don't intersect
    @test !Triangulation.EdgeIntersects(tess, e1011, e1213)
    @test !Triangulation.EdgeIntersects(tess, e1213, e67)

  end

  @testset "FindIntersectingEdges" begin
    local tess = Triangulation.DelaunayTesselation()

    local points = [1.5 1.5; 1.51 1.45; 1.49 1.425; 1.45 1.49; 1.525 1.5; 1.48 1.48; 1.42 1.44]
    push!(tess, points)
    TestHelpers.TestVertexCache(tess)

    #=
    # plot it
    xc, yc = Triangulation.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("findintersecting.svg")
    =#

    # check some intersections for non-existing edges
    intersecting, alreadyIn = Triangulation.FindIntersectingEdges(tess, 6, 8, 1)
    @test !alreadyIn
    @test length(intersecting) == 2
    @test TestHelpers.EdgeConnectsVertices(tess, intersecting[1], 5, 9)
    @test TestHelpers.EdgeConnectsVertices(tess, intersecting[2], 5, 4)

    intersecting, alreadyIn = Triangulation.FindIntersectingEdges(tess, 10, 8, 1)
    @test !alreadyIn
    @test length(intersecting) == 3
    @test TestHelpers.EdgeConnectsVertices(tess, intersecting[1], 6, 9)
    @test TestHelpers.EdgeConnectsVertices(tess, intersecting[2], 5, 9)
    @test TestHelpers.EdgeConnectsVertices(tess, intersecting[3], 5, 4)

    # check some edges that are already in
    intersecting, alreadyIn = Triangulation.FindIntersectingEdges(tess, 6, 9, 1)
    @test alreadyIn
    @test length(intersecting) == 1
    @test TestHelpers.EdgeConnectsVertices(tess, intersecting[1], 6, 9)

    intersecting, alreadyIn = Triangulation.FindIntersectingEdges(tess, 9, 5, 1)
    @test alreadyIn
    @test length(intersecting) == 1
    @test TestHelpers.EdgeConnectsVertices(tess, intersecting[1], 9, 5)
  end

  @testset "RemoveIntersectingEdges" begin
    local tess = Triangulation.DelaunayTesselation()

    local points = [1.5 1.5; 1.51 1.45; 1.49 1.425; 1.45 1.49; 1.525 1.5; 1.48 1.48; 1.42 1.44]
    push!(tess, points)
    TestHelpers.TestVertexCache(tess)

    intersecting, alreadyIn = Triangulation.FindIntersectingEdges(tess, 10, 8, 1)
    @test !alreadyIn
    @test length(intersecting) == 3

    local constraint = Triangulation.RemoveIntersectingEdges!(tess, 10, 8, intersecting)
    TestHelpers.TestVertexCache(tess)
    @test length(intersecting) == 3 # it should not affect the container

    # check that constraint connects the correct vertices and does not
    # intersect any other edges
    @test TestHelpers.EdgeConnectsVertices(tess, constraint, 10, 8)
    for i in 1:4:length(tess.edges)
      if (QuadEdge.Base(i) != QuadEdge.Base(constraint)) && (QuadEdge.Org(tess, i) != QuadEdge.Org(tess, constraint)) && (QuadEdge.Org(tess, i) != QuadEdge.Dest(tess, constraint)) && (QuadEdge.Dest(tess, i) != QuadEdge.Org(tess, constraint)) && (QuadEdge.Dest(tess, i) != QuadEdge.Dest(tess, constraint))
          @test !Triangulation.EdgeIntersects(tess, constraint, i)
      end
    end

    #=
    # plot it
    xc, yc = Triangulation.GetDelaunayCoordinates(tess)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.savefig("removeintersecting.svg")
    =#
  end
end
