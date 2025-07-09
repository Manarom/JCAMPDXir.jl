
using Revise, Plots,OrderedCollections,StaticArrays,DelimitedFiles

include("JCAMPDXir.jl")
test_data_folder = joinpath(".","test","tests data")
python_package_test_data = joinpath(test_data_folder,"jcamp_python")


ASDF_file_decoded = joinpath(test_data_folder,"isopropanol_ASDF_decoded.txt")
asdf_file = joinpath(python_package_test_data,"isopropanol_ASDF.jdx")
asdf_data = readdlm(ASDF_file_decoded)

data = JCAMPDXir.read_jdx_file(asdf_file)
plot(data.x,data.y)
plot!(asdf_data[:,1],asdf_data[:,2])
mat = hcat(asdf_data,data.x,data.y)


files_from_jcamp_py = filter(r->occursin(".jdx",r), readdir(python_package_test_data))
println("\ntesting files from $(python_package_test_data)" )
readed_num = length(files_from_jcamp_py)
problem_files_index = OrderedDict{Int,Pair{Int,String}}()
error_counter = 0

for (i,f) in enumerate(files_from_jcamp_py)
    file_headers = JCAMPDXir.parse_headers(joinpath(python_package_test_data,f))
    try
        @show (i,f)
        data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,f))
    catch ex
        bt = backtrace()
        msg = sprint(showerror, ex, bt)
        global error_counter +=1
        push!(problem_files_index,error_counter=>Pair(i,msg))
    end
end
println("error encountered for $error_counter")
println("error fraction $(100*error_counter/length(files_from_jcamp_py)) %")
problem_files_index[1]
files_from_jcamp_py[28]
data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,files_from_jcamp_py[28])) # reading test file
plot(data.x,data.y)
data.y[end]
data.validation
data.validation.deltax.message


data_headers = JCAMPDXir.parse_headers(joinpath(python_package_test_data,files_from_jcamp_py[6]))

data = JCAMPDXir.read_jdx_file(raw".\test\tests data\jcamp_python\example_multiline_datasets.jdx")
plot(data.x,data.y)
data.headers

T_data = JCAMPDXir.read_jdx_file(raw".\test\tests data\JCAMPDX5.jdx") 
plot(T_data.x,T_data.y)
T_data.headers["NPOINTS"]
length(T_data.y)
head = JCAMPDXir.parse_headers(raw".\test\tests data\JCAMP_XYXY.jdx")
head
include("JCAMPDXir.jl")
jdx_blocks = JCAMPDXir.count_blocks(raw".\test\tests data\jcamp_python\example_compound_file.jdx")
out = JCAMPDXir.read_jdx_file(raw".\test\tests data\jcamp_python\example_compound_file.jdx")
using Plots
plot(out[1].x,out[1].y)
plot!(out[2].x,out[2].y)