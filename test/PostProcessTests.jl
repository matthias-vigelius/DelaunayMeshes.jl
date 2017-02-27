#import Winston

@testset "PostProcessTests" begin
  @testset "refine_valid_triangles" begin
    local mesh = DelaunayMeshes.Mesh()
    DelaunayMeshes.setboundingbox(mesh, [-15.0, 15.0, -15.0, 15.0])
    # seed rng and store seed in file
    local seed = 2628265373
    #local seed = rand(UInt32)
    local fs = open("refine_valid_triangles.seed", "w")
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

    # add constraints
    local innerConstraintVertexList = [x for x in (nvertices+3+nringvertices):-1:(nvertices+3+1)]
    local outerConstraintVertexList = [x for x in (nvertices+3+nringvertices+1):(nvertices + 3 + 2*nringvertices)]
    DelaunayMeshes.addconstraint!(mesh, innerConstraintVertexList)
    DelaunayMeshes.addconstraint!(mesh, outerConstraintVertexList)

    # refine
    DelaunayMeshes.refine_valid_triangles!(mesh)

    # check face locations
    local innerBoundaryEdge = findedgeconnectingvertices(mesh.tesselation, nvertices+3+1, nvertices+3+2)
    @test mesh.tesselation.edges[innerBoundaryEdge].constraintType == DelaunayMeshes.OnLeftSide
    local insideFaces = DelaunayMeshes.getfacesinsideregion(mesh, DelaunayMeshes.rot(innerBoundaryEdge))
    local expLoc = fill(false, length(mesh.tesselation.faces))
    expLoc[insideFaces] = true
    @test expLoc == mesh.faceLocation

    # check conformity of triangles
    check_conformity(mesh.tesselation)

    # check that all triangles are high quality
    for tri in mesh.tesselation
      if (mesh.faceLocation[tri.face])
        @test DelaunayMeshes.check_triangle(mesh, tri)
      end
    end

    # usual checks
    testvertexcache(mesh.tesselation)

    #=
    # plot it
    xc, yc = DelaunayMeshes.getdelaunaycoordinates(mesh.tesselation)
    p = Winston.FramedPlot(aspect_ratio=1
      #,xrange=[1.496 1.50], yrange=[1.256, 1.26]
      )
    Winston.add(p, Winston.Curve(xc, yc))
    Winston.savefig(p, "refine_valid_triangles.svg")
    =#


  end

  @testset "refine_triangle" begin
    local mesh = DelaunayMeshes.Mesh()
    DelaunayMeshes.setboundingbox(mesh, [-15.0, 15.0, -15.0, 15.0])

    push!(mesh, [-5.0 -5.0; 5.0 -5.0; 0.0 5.0])

    # find face foce for vertices
    local testTriangle = DelaunayMeshes.Triangle(mesh.tesselation, 4, 5, 6)

    # refine
    DelaunayMeshes.refine_triangle(mesh, testTriangle)

    # check that there is a new vertex and that it's at the offcenter
    @test length(mesh.tesselation.vertices) == 7
    local offcenter = DelaunayMeshes.compute_offcenter(mesh, 4, 5, 6)
    @test_approx_eq(mesh.tesselation.vertices[7]', offcenter')

    #=
    # plot it
    xc, yc = getplotedges(mesh.tesselation)
    local vertexPoints = reduce(hcat, mesh.tesselation.vertices)'
    p = Winston.FramedPlot(aspect_ratio=1
      ,xrange=[1.45, 1.55], yrange=[1.2, 1.3]
      )
    Winston.add(p, Winston.Points(vertexPoints[:,1], vertexPoints[:,2], kind="circle", color="red"))
    Winston.add(p, Winston.Curve(xc, yc))
    Winston.savefig(p, "refine_triangle.svg")
    =#
  end

  @testset "refine_triangle_constraints" begin
    local mesh = DelaunayMeshes.Mesh()
    DelaunayMeshes.setboundingbox(mesh, [-15.0, 15.0, -15.0, 15.0])

    push!(mesh, [-5.0 -5.0; 5.0 -5.0; 0.0 5.0])

    # find face foce for vertices
    local tri = DelaunayMeshes.Triangle(mesh.tesselation, 4, 5, 6)

    # set edges as constraints
    DelaunayMeshes.setconstraint!(mesh.tesselation, tri.edges[1], DelaunayMeshes.OnLeftSide)
    DelaunayMeshes.setconstraint!(mesh.tesselation, tri.edges[2], DelaunayMeshes.OnLeftSide)

    # refine
    DelaunayMeshes.refine_triangle(mesh, tri)

    # check that there is a new vertex and that it's at the edge centers
    @test length(mesh.tesselation.vertices) == 8
    midpoint1 = 0.5*(mesh.tesselation.vertices[4] + mesh.tesselation.vertices[5])
    midpoint2 = 0.5*(mesh.tesselation.vertices[5] + mesh.tesselation.vertices[6])
    @test_approx_eq(mesh.tesselation.vertices[7]', midpoint1')
    @test_approx_eq(mesh.tesselation.vertices[8]', midpoint2')

    # check that the edges going to and from these vertices are marked as
    # segments
    local e47 = findedgeconnectingvertices(mesh.tesselation, 4, 7)
    local e75 = findedgeconnectingvertices(mesh.tesselation, 7, 5)
    @test mesh.tesselation.edges[e47].constraintType == DelaunayMeshes.OnLeftSide
    @test mesh.tesselation.edges[e75].constraintType == DelaunayMeshes.OnLeftSide
    @test mesh.tesselation.edges[DelaunayMeshes.sym(e47)].constraintType == DelaunayMeshes.OnRightSide
    @test mesh.tesselation.edges[DelaunayMeshes.sym(e75)].constraintType == DelaunayMeshes.OnRightSide
    local e58 = findedgeconnectingvertices(mesh.tesselation, 5, 8)
    local e86 = findedgeconnectingvertices(mesh.tesselation, 8, 6)
    @test mesh.tesselation.edges[e58].constraintType == DelaunayMeshes.OnLeftSide
    @test mesh.tesselation.edges[e86].constraintType == DelaunayMeshes.OnLeftSide
    @test mesh.tesselation.edges[DelaunayMeshes.sym(e58)].constraintType == DelaunayMeshes.OnRightSide
    @test mesh.tesselation.edges[DelaunayMeshes.sym(e86)].constraintType == DelaunayMeshes.OnRightSide

    #=
    # plot it
    xc, yc = getplotedges(mesh.tesselation)
    local vertexPoints = reduce(hcat, mesh.tesselation.vertices)'
    p = Winston.FramedPlot(aspect_ratio=1
      ,xrange=[1.45, 1.55], yrange=[1.2, 1.3]
      )
    Winston.add(p, Winston.Points(vertexPoints[:,1], vertexPoints[:,2], kind="circle", color="red"))
    Winston.add(p, Winston.Curve(xc, yc))
    Winston.savefig(p, "refine_triangle.svg")
    =#
  end


  @testset "get_shortest_edge" for par in [[-1.0, 2.0, 1.3, 5.0, 7.0, 7.5], [1.0, -2.0, 3.141, 10.0, 12.0, 21.5], [0.0, 0.0, 0.0, 1.0, 2.0, 2.5]] begin
    cx, cy, alpha0, a, b, c = par

    mesh = DelaunayMeshes.Mesh()

    s = 0.5*(a+b+c)
    r = sqrt((s-a)*(s-b)*(s-c)/s)
    alpha = 2.0*(atan(r/(s-a)))
    beta = 2.0*(atan(r/(s-b)))
    gamma = 2.0*(atan(r/(s-c)))
    (v1x, v1y) = (cx + b * cos(alpha + alpha0), cy + b * sin(alpha + alpha0))
    (v2x, v2y) = (cx + c * cos(alpha0), cy + c * sin(alpha0))
    mesh.tesselation.vertices = [[cx, cy],[v2x, v2y],[v1x, v1y]]

    p, q, r = DelaunayMeshes.get_shortest_edge(mesh, 1, 2, 3)
    @test p == 1
    @test q == 2
    @test r == 3
    p, q, r = DelaunayMeshes.get_shortest_edge(mesh, 2, 3, 1)
    @test p == 1
    @test q == 2
    @test r == 3
    p, q, r = DelaunayMeshes.get_shortest_edge(mesh, 3, 1, 2)
    @test p == 1
    @test q == 2
    @test r == 3
  end
  end

  @testset "compute_circumcenter_radius" for r in [0.1, 1.1], c in [[1.6, 1.3], [1.25, 1.7]], t in [[0, pi/2., pi], [pi/8, pi-pi/16, pi+pi/16]] begin
    # generate triangle on given circumcircle
    local cx = c[1]
    local cy = c[2]
    local t1 = t[1]
    local t2 = t[2]
    local t3 = t[3]
    mesh = create_mesh_with_triangle(cx, cy, r, t1, t2, t3)

    sx, sy, rs = DelaunayMeshes.compute_circumcenter_radius(mesh, 1, 2, 3)

    @test(abs(sx-cx) <= 1e-4)
    @test(abs(sy-cy) <= 1e-4)
    @test(abs(rs-r*r) <= 1e-4)
  end
  end

  @testset "compute_offcenter" for c in [[1.6, 1.3], [1.25, 1.7]] , a0 in [pi/24., pi/16., pi/8.], offset in [0., pi/3., pi, 1.5*pi] begin
    local cx = c[1]
    local cy = c[2]
    local r = 2.
    local t1 = pi - a0 + offset
    local t2 = pi + a0 + offset
    local t3 = pi/4. + offset
    mesh = create_mesh_with_triangle(cx, cy, r, t1, t2, t3)

    local offcenter = DelaunayMeshes.compute_offcenter(mesh, 1, 2, 3)

    # check that new triangle p - q - oc has correct edge/circumcircle radius
    mesh.tesselation.vertices[3] = offcenter
    sx, sy, rs = DelaunayMeshes.compute_circumcenter_radius(mesh, 1, 2, 3)
    local lpq = norm(mesh.tesselation.vertices[1] - mesh.tesselation.vertices[2])
    @test sqrt(rs)/lpq <= DelaunayMeshes.beta
  end
  end

  @testset "vertex_encroaches_segment" for rd in [0.8, 1.2] begin
    # seed rng and store seed in file
    #local seed =2410381395
    local seed = rand(UInt32)
    local fs = open("vertex_encroaches_segment.seed", "w")
    write(fs, string(seed))
    close(fs)
    srand(seed)

    mesh = DelaunayMeshes.Mesh()

    # mark first edge as boundary
    DelaunayMeshes.setconstraint!(mesh.tesselation, 1, DelaunayMeshes.OnLeftSide)
    mesh.tesselation.edges[5].constraintType = DelaunayMeshes.NoConstraint
    mesh.tesselation.edges[7].constraintType = DelaunayMeshes.NoConstraint

    # create potentially encroaching vertex for that edge
    me = 0.5*(mesh.tesselation.vertices[1] + mesh.tesselation.vertices[2])
    l = norm(mesh.tesselation.vertices[1] - mesh.tesselation.vertices[2])
    θ = rand(Float64) * 2. * pi
    vp = convert(DelaunayMeshes.Point, me + rd * 0.5 * l * [cos(θ), sin(θ)])
    @test DelaunayMeshes.vertex_encroaches_segment(mesh, 1, vp) == (rd < 1.0)
    @test DelaunayMeshes.vertex_encroaches_segment(mesh, 3, vp) == (rd < 1.0)

    # create possibly encroaching vertex for non boundary
    me = 0.5*(mesh.tesselation.vertices[2] + mesh.tesselation.vertices[3])
    l = norm(mesh.tesselation.vertices[2] - mesh.tesselation.vertices[3])
    θ = rand(Float64) * 2. * pi
    vp = convert(DelaunayMeshes.Point, me + rd * 0.5 * l * [cos(θ), sin(θ)])
    @test DelaunayMeshes.vertex_encroaches_segment(mesh, 5, vp) == false
    @test DelaunayMeshes.vertex_encroaches_segment(mesh, 7, vp) == false
  end
  end
end
