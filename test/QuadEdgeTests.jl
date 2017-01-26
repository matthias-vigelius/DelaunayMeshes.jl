@testset "QuadEdgeTests" begin
  @testset "base" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test base(baseIndex) == baseIndex
    @test base(baseIndex+1) == baseIndex
    @test base(baseIndex+2) == baseIndex
    @test base(baseIndex+3) == baseIndex
  end

  @testset "sym" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test sym(baseIndex) == baseIndex + 2
    @test sym(baseIndex+2) == baseIndex
    @test sym(baseIndex+1) == baseIndex + 3
    @test sym(baseIndex+3) == baseIndex + 1
  end

  @testset "rot" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test rot(baseIndex) == baseIndex + 1
    @test rot(baseIndex+1) == baseIndex + 2
    @test rot(baseIndex+2) == baseIndex + 3
    @test rot(baseIndex+3) == baseIndex
  end

  @testset "invrot" begin
    local baseIndex = rand(0:20) * 4 + 1
    @test invrot(rot(baseIndex)) == baseIndex
    @test invrot(rot(baseIndex+1)) == baseIndex + 1
    @test invrot(rot(baseIndex+2)) == baseIndex + 2
    @test invrot(rot(baseIndex+3)) == baseIndex + 3
  end

  @testset "makeedge!" begin
    # make empty subdivision and check if make edge adds edges
    local sd = SubDivision{Int}([],[],[],[])
    @test length(sd.edges)==0

    local ne = makeedge!(sd)
    @test ne == 1

    @test sd.edges[1].origin == 0
    @test sd.edges[1].next == 1
    @test sd.edges[1].constraintType == NoConstraint
    @test sd.edges[2].origin == 0
    @test sd.edges[2].next == 4
    @test sd.edges[2].constraintType == NoConstraint
    @test sd.edges[3].origin == 0
    @test sd.edges[3].next == 3
    @test sd.edges[3].constraintType == NoConstraint
    @test sd.edges[4].origin == 0
    @test sd.edges[4].next == 2
    @test sd.edges[4].constraintType == NoConstraint
  end

  @testset "splice!" begin
    local sd = SubDivision{Int}([],[],[],[])
    local ne1 = makeedge!(sd)
    local ne2 = makeedge!(sd)

    # splice edges together to form a single ring
    splice!(sd, ne1, ne2)
    @test onext(sd, ne1) == ne2
    @test onext(sd, ne2) == ne1
    @test onext(sd, ne1 + 1) == ne2 + 3
    @test onext(sd, ne1 + 3) == ne1 + 1
    @test onext(sd, ne2 + 1) == ne1 + 3
    @test onext(sd, ne2 + 3) == ne2 + 1

    # cut rings again
    splice!(sd, ne1, ne2)
    @test onext(sd, ne1) == ne1
    @test onext(sd, ne2) == ne2
    @test onext(sd, ne1 + 1) == ne1 + 3
    @test onext(sd, ne1 + 3) == ne1 + 1
    @test onext(sd, ne2 + 1) == ne2 + 3
    @test onext(sd, ne2 + 3) == ne2 + 1
  end


  @testset "endpoints" begin
    local sd = SubDivision{Float64}([],[],[],[])
    makeedge!(sd)
    sd.vertexCache = [-1,-1,-1]

    # check setting of end points works correctly
    endpoints!(sd, 1,1,2,3,4)
    @test org(sd, 1) == 1
    @test dest(sd, 1) == 2
    @test org(sd, sym(1)) == 2
    @test dest(sd, sym(1)) == 1
    @test sd.vertexCache[1] == 1
    @test sd.vertexCache[2] == 3
    @test lface(sd, 1) == 3
    @test rface(sd, 1) == 4
    @test lface(sd, sym(1)) == 4
    @test rface(sd, sym(1)) == 3

    # check end points works correctly for sym edge
    endpoints!(sd, 3,1,2,3,4)
    @test sd.edges[1].origin == 2
    @test sd.edges[3].origin == 1
    @test sd.vertexCache[1] == 3
    @test sd.vertexCache[2] == 1
    @test lface(sd, 3) == 3
    @test rface(sd, 3) == 4
    @test lface(sd, sym(3)) == 4
    @test rface(sd, sym(3)) == 3
  end

  @testset "subdivision" begin
    local sd = createsubdivision(7.0,3.0,5.0, 1.0, 2.0)
    @test isa(sd, SubDivision{Float64})
    TestHelpers.TestVertexCache(sd)

    # test that created triangle works
    local ea = 1
    @test org(sd, ea) == 1
    @test dest(sd, ea) == 2
    @test org(sd, rot(ea)) == 2
    @test dest(sd, rot(ea)) == 1

    local eb = lnext(sd, ea)
    @test org(sd, eb) == 2
    @test dest(sd, eb) == 3
    @test org(sd, rot(eb)) == 2
    @test dest(sd, rot(eb)) == 1

    local ec = lnext(sd, eb)
    @test org(sd, ec) == 3
    @test dest(sd, ec) == 1
    @test org(sd, rot(ec)) == 2
    @test dest(sd, rot(ec)) == 1
  end

  @testset "connect" begin
    local sd = createsubdivision(7.0,3.0,5.0, 1.0, 2.0)
    local ea = 1
    local eb = 5
    local ec = 9
    local innerFaceIndex = 1
    local outerFaceIndex = 2
    local newFaceIndex = 3

    # add single edge to main triangle starting at vertex 1
    push!(sd.vertices, 8.0)
    push!(sd.vertexCache, -1.)
    local ed = makeedge!(sd)
    endpoints!(sd, ed, 1, 4, innerFaceIndex, innerFaceIndex)
    splice!(sd, ed, ea)

    # and connect it to vertex 2
    local ee = connect!(sd, ea, sym(ed), 3.0)

    @test org(sd, rot(ea)) == outerFaceIndex
    @test dest(sd, rot(ea)) == innerFaceIndex
    @test org(sd, rot(eb)) == outerFaceIndex
    @test dest(sd, rot(eb)) == newFaceIndex
    @test org(sd, rot(ec)) == outerFaceIndex
    @test dest(sd, rot(ec)) == newFaceIndex
    @test org(sd, rot(ed)) == innerFaceIndex
    @test dest(sd, rot(ed)) == newFaceIndex
    @test org(sd, rot(ee)) == newFaceIndex
    @test dest(sd, rot(ee)) == innerFaceIndex
  end
end
