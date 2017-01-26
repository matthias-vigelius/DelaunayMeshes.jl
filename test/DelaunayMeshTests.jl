using Base.Test

import DelaunayMesh
import TestHelpers
import DiffEqPDEBase

import Winston

@testset "DelaunayMeshTests" begin
  @testset "AddConstraint" begin
    local mesh = DelaunayMesh.Mesh()
    DelaunayMesh.SetBoundingBox(mesh, [-15.0, 15.0, -15.0, 15.0])
    # seed rng and store seed in file
    local seed = 4139605370
    #local seed = rand(UInt32)
    local fs = open("AddConstraint.seed", "w")
    write(fs, string(seed))
    close(fs)

    # push some random vertices
    srand(seed)
    local nvertices = 200
    local points = rand(Float64, nvertices, 2)*30. - 15.
    push!(mesh, points)

    local nringvertices = 40

    # push inner ring constraint vertices
    local θ = linspace(0., 2.*π, nringvertices + 1)[1:end-1]
    local rinner = 5.
    local router = 10.
    local innerRing = rinner * [cos(θ'); sin(θ')]'
    local outerRing = router * [cos(θ'); sin(θ')]'

    push!(mesh, innerRing)
    push!(mesh, outerRing)

    local faultyTriangles = TestHelpers.testdelaunayness(mesh.tesselation)
    @test length(faultyTriangles) == 0
    @test length(mesh.tesselation.vertices) == (nvertices + 3 + 2*nringvertices)

    # add constraints
    local innerConstraintVertexList = [x for x in (nvertices+3+nringvertices):-1:(nvertices+3+1)]
    local outerConstraintVertexList = [x for x in (nvertices+3+nringvertices+1):(nvertices + 3 + 2*nringvertices)]

    DelaunayMesh.AddConstraint!(mesh, innerConstraintVertexList)
    DelaunayMesh.AddConstraint!(mesh, outerConstraintVertexList)

    # check that constraints are in
    local innerConstrMatrix = [innerConstraintVertexList'; circshift(innerConstraintVertexList, -1)']'
    local outerConstrMatrix = [outerConstraintVertexList'; circshift(outerConstraintVertexList, -1)']'
    local constraints = [innerConstrMatrix; outerConstrMatrix]

    TestHelpers.validateconstraints(mesh.tesselation, constraints)
    TestHelpers.testvertexcache(mesh.tesselation)

    # get voronoiVertices
    local vorVert = DelaunayMesh.GetVoronoiVertices(mesh)

    # mark all that are between inner and outer circle
    local unScaledVorVert = DelaunayMesh.UnScalePoints(mesh, vorVert)
    local unScaledSquare = unScaledVorVert.*unScaledVorVert
    local rVor = sqrt(unScaledSquare[:,1] + unScaledSquare[:,2])
    local interiorVerticesBool = (rVor .>= rinner * ones(rVor)) & (rVor .<= router * ones(rVor))

    # plot it
    xc, yc = Triangulation.GetDelaunayCoordinates(mesh.tesselation)
    p = Winston.FramedPlot(aspect_ratio=1
      #,xrange=[1.496 1.50], yrange=[1.256, 1.26]
      )
    Winston.add(p, Winston.Curve(xc, yc))
    #Winston.add(p, Winston.Points(vorVert[interiorVerticesBool,1], vorVert[interiorVerticesBool,2], kind="circle", color="green"))
    Winston.add(p, Winston.Points(
      vorVert[mesh.faceLocation,1], vorVert[mesh.faceLocation,2],
      kind="circle", color="red"))
    Winston.savefig(p, "AddConstraint.svg")

    # check that this agrees with the interior array of the triangulation
    @test length(interiorVerticesBool) == length(mesh.faceLocation)
    if !(interiorVerticesBool == mesh.faceLocation)
      for (i,e) in enumerate(interiorVerticesBool)
        if (e != mesh.faceLocation[i])
          local pos = unScaledVorVert[i, :]
          local r = sqrt(pos[1]^2 + pos[2]^2)
          println("Pos $i: expecting $e, got $(mesh.faceLocation[i]). $pos and $r (scaled:$(vorVert[i,:]))")
        end
      end
    end

    @test interiorVerticesBool == mesh.faceLocation
  end

  @testset "FindConstrainedTriangles" begin
    local mesh = DelaunayMesh.Mesh()

    # seed rng and store seed in file
    #local seed = 2572853444
    local seed = rand(UInt32)
    local fs = open("FindConstrainedTriangles.seed", "w")
    write(fs, string(seed))
    close(fs)

    # push some random vertices
    srand(seed)
    local nvertices = 20
    local points = rand(Float64, nvertices, 2)*5.0 - 2.5
    push!(mesh, points)

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
    push!(mesh, constraintVertices)

    local faultyTriangles = TestHelpers.testdelaunayness(mesh.tesselation)
    #=
    if (length(faultyTriangles) > 0)
      for f in faultyTriangles
        local ftIndex = [f[1], f[2], f[3], f[4]]
        local vertexPoints = reduce(hcat, mesh.tesselation.vertices[ftIndex])'
        PyPlot.plot(vertexPoints[:,1], vertexPoints[:,2], "o")
      end
    end
    =#
    @test length(faultyTriangles) == 0
    @test length(mesh.tesselation.vertices) == (nvertices + 8 + 3)

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
      Triangulation.DelaunayMeshes.insertconstraint!(mesh.tesselation, constraints[i, 1], constraints[i,2])
    end

    TestHelpers.validateconstraints(mesh.tesselation, constraints)
    TestHelpers.testvertexcache(mesh.tesselation)

    # get voronoiVertices
    local vorVert = DelaunayMesh.GetVoronoiVertices(mesh)

    # get faces in connected region
    local startEdge = TestHelpers.findedgeconnectingvertices(mesh.tesselation, n+1, n+2)
    local startDualEdge = sym(Rot(startEdge))
    local connectedVertices = DelaunayMesh.GetFacesInsideRegion(mesh, startDualEdge)
    local connectedBool = fill(false, size(mesh.tesselation.faces, 1))
    connectedBool[connectedVertices] = true

    # get all faces inside constraint for comparison
    local insideConstraint = TestHelpers.PointsInConstrainedRegion(mesh, vorVert)

    @test connectedBool == insideConstraint

    # plot it
    #=
    xc, yc = Triangulation.GetDelaunayCoordinates(mesh.tesselation)
    PyPlot.clf()
    PyPlot.plot(xc, yc)
    PyPlot.plot(get(mesh.voronoiVertices)[insideConstraint,1], get(mesh.voronoiVertices)[insideConstraint,2], "o")
    PyPlot.savefig("FindConstrainedTriangles.svg")
    =#

  end

  @testset "push" begin
    local mesh = DelaunayMesh.Mesh()

    # seed rng and store seed in file
    #local seed = 1985005665
    local seed = rand(UInt32)
    local fs = open("DelaunayMeshPushRandom.seed", "w")
    write(fs, string(seed))
    close(fs)

    # get dimensions of bounding box
    local testBB = [rand(), rand(), rand(), rand()] * 1e5 - 0.5e4

    # push some random vertices
    srand(seed)
    local nvertices = 200
    local points = TestHelpers.GetRandomNumbersInsideBoundingBox(testBB, nvertices)
    push!(mesh, points)

    # check they are inside bounding box
    local vertexPoints = reduce(hcat, mesh.tesselation.vertices[4:end])'
    local minX = minimum(vertexPoints[:,1])
    local maxX = maximum(vertexPoints[:,1])
    local minY = minimum(vertexPoints[:,2])
    local maxY = maximum(vertexPoints[:,2])

    local dh = 0.05/2.
    @test minX >= (1.5 - dh)
    @test maxX <= (1.5 + dh)
    @test minY >= (1.25 - dh)
    @test maxY <= (1.25 + dh)

    # push some more vertices that are guaranteed to be inside
    local smallBB = [
      minimum(points[:,1]), maximum(points[:,1]),
      minimum(points[:,2]), maximum(points[:,2]),
    ]

    local morePoints = TestHelpers.GetRandomNumbersInsideBoundingBox(smallBB, nvertices)
    push!(mesh, morePoints)

    # check they are still inside bounding box
    vertexPoints = reduce(hcat, mesh.tesselation.vertices[4:end])'
    minX = minimum(vertexPoints[:,1])
    maxX = maximum(vertexPoints[:,1])
    minY = minimum(vertexPoints[:,2])
    maxY = maximum(vertexPoints[:,2])

    @test minX >= (1.5 - dh)
    @test maxX <= (1.5 + dh)
    @test minY >= (1.25 - dh)
    @test maxY <= (1.25 + dh)

    # check that the vertex locations are correct
    @test length(mesh.faceLocation) == length(mesh.tesselation.faces)
    @test mesh.faceLocation[1]
    @test !mesh.faceLocation[2]
    @test all(mesh.faceLocation[3:end])
  end
end
