using JCAMPDXir
using Test
test_data_folder = joinpath(@__DIR__(),"tests data")
fname = joinpath(test_data_folder,"JCAMP_test_file.jdx")
@testset "JCAMPDXir.jl" begin
    # Write your tests here.
    data = JCAMPDXir.read_jdx_file(fname)
end
