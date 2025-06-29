
#test_data_folder = joinpath(".","test","tests data")
#python_package_test_data = joinpath(test_data_folder,"jcamp_python")
include("JCAMPDXir.jl")
#data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,"1,2-dichloroethane.jdx"))
global jdx_blocks
jdx_blocks = JCAMPDXir.count_blocks(raw".\test\tests data\jcamp_python\example_compound_file.jdx")