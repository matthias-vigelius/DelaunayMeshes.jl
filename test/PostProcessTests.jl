#TODO: check also case where offcenter is to the right of circum center (a0 = 16. in asy..)
@testset "PostProcessTests" begin
  @testset "compute_circumcenter_radius" for r in [0.1, 1.1], c in [[1.6, 1.3], [1.25, 1.7]], t in [[0, pi/2., pi], [pi/8, pi-pi/16, pi+pi/16]] begin
    # generate triangle on given circumcircle
    local cx = c[1]
    local cy = c[2]
    local t1 = t[1]
    local t2 = t[2]
    local t3 = t[3]
    local v1 = [cx + r * cos(t1), cy + r * sin(t1)]
    local v2 = [cx + r * cos(t2), cy + r * sin(t2)]
    local v3 = [cx + r * cos(t3), cy + r * sin(t3)]
    local mesh = DelaunayMeshes.Mesh()
    mesh.tesselation.vertices = [v1, v2, v3]
    sx, sy, rs = DelaunayMeshes.compute_circumcenter_radius(mesh, 1, 2, 3)

    @test(abs(sx-cx) <= 1e-4)
    @test(abs(sy-cy) <= 1e-4)
    @test(abs(rs-r*r) <= 1e-4)
  end
  end
end
