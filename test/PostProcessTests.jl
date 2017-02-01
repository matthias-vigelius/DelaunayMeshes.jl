#TODO: check also case where offcenter is to the right of circum center (a0 = 16. in asy..)
@testset "PostProcessTests" begin
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
