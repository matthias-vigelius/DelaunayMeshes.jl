@testset "QuadEdgeTests" begin
  @testset "base" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test DelaunayMeshes.base(baseIndex) == baseIndex
    @test DelaunayMeshes.base(baseIndex+1) == baseIndex
    @test DelaunayMeshes.base(baseIndex+2) == baseIndex
    @test DelaunayMeshes.base(baseIndex+3) == baseIndex
  end

  @testset "sym" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test DelaunayMeshes.sym(baseIndex) == baseIndex + 2
    @test DelaunayMeshes.sym(baseIndex+2) == baseIndex
    @test DelaunayMeshes.sym(baseIndex+1) == baseIndex + 3
    @test DelaunayMeshes.sym(baseIndex+3) == baseIndex + 1
  end

  @testset "rot" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test DelaunayMeshes.rot(baseIndex) == baseIndex + 1
    @test DelaunayMeshes.rot(baseIndex+1) == baseIndex + 2
    @test DelaunayMeshes.rot(baseIndex+2) == baseIndex + 3
    @test DelaunayMeshes.rot(baseIndex+3) == baseIndex
  end

  @testset "invrot" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test DelaunayMeshes.invrot(DelaunayMeshes.rot(baseIndex)) == baseIndex
    @test DelaunayMeshes.invrot(DelaunayMeshes.rot(baseIndex+1)) == baseIndex + 1
    @test DelaunayMeshes.invrot(DelaunayMeshes.rot(baseIndex+2)) == baseIndex + 2
    @test DelaunayMeshes.invrot(DelaunayMeshes.rot(baseIndex+3)) == baseIndex + 3
  end

  @testset "makeedge!" begin
    # make empty subdivision and check if make edge adds edges
    local sd = DelaunayMeshes.SubDivision{Int}([],[],[],[])
    @test length(sd.edges)==0

    local ne = DelaunayMeshes.makeedge!(sd)
    @test ne == 1

    @test sd.edges[1].origin == 0
    @test sd.edges[1].next == 1
    @test sd.edges[1].constraintType == DelaunayMeshes.NoConstraint
    @test sd.edges[2].origin == 0
    @test sd.edges[2].next == 4
    @test sd.edges[2].constraintType == DelaunayMeshes.NoConstraint
    @test sd.edges[3].origin == 0
    @test sd.edges[3].next == 3
    @test sd.edges[3].constraintType == DelaunayMeshes.NoConstraint
    @test sd.edges[4].origin == 0
    @test sd.edges[4].next == 2
    @test sd.edges[4].constraintType == DelaunayMeshes.NoConstraint
  end

  @testset "splice!" begin
    local sd = DelaunayMeshes.SubDivision{Int}([],[],[],[])
    local ne1 = DelaunayMeshes.makeedge!(sd)
    local ne2 = DelaunayMeshes.makeedge!(sd)

    # splice edges together to form a single ring
    DelaunayMeshes.splice!(sd, ne1, ne2)
    @test DelaunayMeshes.onext(sd, ne1) == ne2
    @test DelaunayMeshes.onext(sd, ne2) == ne1
    @test DelaunayMeshes.onext(sd, ne1 + 1) == ne2 + 3
    @test DelaunayMeshes.onext(sd, ne1 + 3) == ne1 + 1
    @test DelaunayMeshes.onext(sd, ne2 + 1) == ne1 + 3
    @test DelaunayMeshes.onext(sd, ne2 + 3) == ne2 + 1

    # cut rings again
    DelaunayMeshes.splice!(sd, ne1, ne2)
    @test DelaunayMeshes.onext(sd, ne1) == ne1
    @test DelaunayMeshes.onext(sd, ne2) == ne2
    @test DelaunayMeshes.onext(sd, ne1 + 1) == ne1 + 3
    @test DelaunayMeshes.onext(sd, ne1 + 3) == ne1 + 1
    @test DelaunayMeshes.onext(sd, ne2 + 1) == ne2 + 3
    @test DelaunayMeshes.onext(sd, ne2 + 3) == ne2 + 1
  end


  @testset "endpoints" begin
    local sd = DelaunayMeshes.SubDivision{Float64}([],[],[],[])
    DelaunayMeshes.makeedge!(sd)
    sd.vertexCache = [-1,-1,-1]

    # check setting of end points works correctly
    DelaunayMeshes.endpoints!(sd, 1,1,2,3,4)
    @test DelaunayMeshes.org(sd, 1) == 1
    @test DelaunayMeshes.dest(sd, 1) == 2
    @test DelaunayMeshes.org(sd, DelaunayMeshes.sym(1)) == 2
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.sym(1)) == 1
    @test sd.vertexCache[1] == 1
    @test sd.vertexCache[2] == 3
    @test DelaunayMeshes.lface(sd, 1) == 3
    @test DelaunayMeshes.rface(sd, 1) == 4
    @test DelaunayMeshes.lface(sd, DelaunayMeshes.sym(1)) == 4
    @test DelaunayMeshes.rface(sd, DelaunayMeshes.sym(1)) == 3

    # check end points works correctly for sym edge
    DelaunayMeshes.endpoints!(sd, 3,1,2,3,4)
    @test sd.edges[1].origin == 2
    @test sd.edges[3].origin == 1
    @test sd.vertexCache[1] == 3
    @test sd.vertexCache[2] == 1
    @test DelaunayMeshes.lface(sd, 3) == 3
    @test DelaunayMeshes.rface(sd, 3) == 4
    @test DelaunayMeshes.lface(sd, DelaunayMeshes.sym(3)) == 4
    @test DelaunayMeshes.rface(sd, DelaunayMeshes.sym(3)) == 3
  end

  @testset "subdivision" begin
    local sd = DelaunayMeshes.create_subdivision(7.0,3.0,5.0, 1.0, 2.0)
    @test isa(sd, DelaunayMeshes.SubDivision{Float64})
    testvertexcache(sd)

    # test that created triangle works
    local ea = 1
    @test DelaunayMeshes.org(sd, ea) == 1
    @test DelaunayMeshes.dest(sd, ea) == 2
    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(ea)) == 2
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(ea)) == 1

    local eb = DelaunayMeshes.lnext(sd, ea)
    @test DelaunayMeshes.org(sd, eb) == 2
    @test DelaunayMeshes.dest(sd, eb) == 3
    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(eb)) == 2
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(eb)) == 1

    local ec = DelaunayMeshes.lnext(sd, eb)
    @test DelaunayMeshes.org(sd, ec) == 3
    @test DelaunayMeshes.dest(sd, ec) == 1
    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(ec)) == 2
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(ec)) == 1
  end

  @testset "connect" begin
    local sd = DelaunayMeshes.create_subdivision(7.0,3.0,5.0, 1.0, 2.0)
    local ea = 1
    local eb = 5
    local ec = 9
    local innerFaceIndex = 1
    local outerFaceIndex = 2
    local newFaceIndex = 3

    # add single edge to main triangle starting at vertex 1
    push!(sd.vertices, 8.0)
    push!(sd.vertexCache, -1.)
    local ed = DelaunayMeshes.makeedge!(sd)
    DelaunayMeshes.endpoints!(sd, ed, 1, 4, innerFaceIndex, innerFaceIndex)
    DelaunayMeshes.splice!(sd, ed, ea)

    # and connect it to vertex 2
    local ee = DelaunayMeshes.connect!(sd, ea, DelaunayMeshes.sym(ed), 3.0)

    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(ea)) == outerFaceIndex
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(ea)) == innerFaceIndex
    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(eb)) == outerFaceIndex
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(eb)) == newFaceIndex
    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(ec)) == outerFaceIndex
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(ec)) == newFaceIndex
    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(ed)) == innerFaceIndex
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(ed)) == newFaceIndex
    @test DelaunayMeshes.org(sd, DelaunayMeshes.rot(ee)) == newFaceIndex
    @test DelaunayMeshes.dest(sd, DelaunayMeshes.rot(ee)) == innerFaceIndex
  end
end
