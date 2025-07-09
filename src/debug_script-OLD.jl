
#test_data_folder = joinpath(".","test","tests data")
#python_package_test_data = joinpath(test_data_folder,"jcamp_python")
include("JCAMPDXir.jl")
#data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,"1,2-dichloroethane.jdx"))
test_data_folder = joinpath(".","test","tests data")
python_package_test_data = joinpath(test_data_folder,"jcamp_python")
ASDF_file_decoded = joinpath(test_data_folder,"isopropanol_ASDF_decoded.txt")
asd_file = joinpath(python_package_test_data,"isopropanol_ASDF.jdx")
#asd_file = joinpath(python_package_test_data,"ethanol2.jdx")

data = JCAMPDXir.read_jdx_file(asd_file)
using Plots
plot(data.x,data.y)

data.validation
data.validation