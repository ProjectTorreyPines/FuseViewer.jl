using FuseViewer
using Test

@testset "FuseViewer" begin
    for name in (:launch, :open_dir, :start, :wait_for_close)
        @test isdefined(FuseViewer, name)
        @test Base.isexported(FuseViewer, name)
    end
    @test_throws Exception FuseViewer.open_dir(mktempdir(); display_window=false)
end
