#debug script
using Revise, Plots,OrderedCollections,StaticArrays
test_data_folder = joinpath(".","test","tests data")
python_package_test_data = joinpath(test_data_folder,"jcamp_python")
includet("JCAMPDXir.jl")

line = "575.17-3042244 1597332-970474 1254921-1092859 2206136"
jdx = JCAMPDXir.JDXblock()
JCAMPDXir.addline!(JCAMPDXir.XYYline,jdx,line)
jdx.y_data
data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,"benzene.jdx")) # reading test file
plot(data.x,data.y)
a = Float64[]
JCAMPDXir.split_PAC_string!(a,1,"+345-23.0-4")
a
s = MVector{3}(a)
JCAMPDXir.split_PAC_string!(a,3,"+345-23.0-4")
s


line = "11995.21,    32112;   11991.36,    32505;   11987.5,    32727;   11983.64,    33481;"
jdx = JCAMPDXir.JDXblock()
JCAMPDXir.addline!(JCAMPDXir.XYXYline,jdx,line)
jdx.x_data
jdx.y_data
files_from_jcamp_py = filter(r->occursin(".jdx",r), readdir(python_package_test_data))
println("\ntesting files from $(python_package_test_data)" )
readed_num = length(files_from_jcamp_py)

problem_files_index = OrderedDict{Int,Pair{Int,String}}()
error_counter = 0

for (i,f) in enumerate(files_from_jcamp_py)
    file_headers = JCAMPDXir.parse_headers(joinpath(python_package_test_data,f))
    if haskey(file_headers,"JCAMP-DX") && isa(file_headers["JCAMP-DX"],Number) && file_headers["JCAMP-DX"]<=4.24
        @show file_headers["JCAMP-DX"]
    else
        continue
    end
    try
        data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,f))
    catch ex
        bt = backtrace()
        msg = sprint(showerror, ex, bt)
        global error_counter +=1
        push!(problem_files_index,error_counter=>Pair(i,msg))
    end
end
problem_files_index[1]
error_counter
println("error $(100*error_counter/length(files_from_jcamp_py)) %")
files_from_jcamp_py[43]
data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,files_from_jcamp_py[43])) # reading test file
plot(data.x,data.y)
data.y[end]

data_headers = JCAMPDXir.parse_headers(joinpath(python_package_test_data,files_from_jcamp_py[6]))

data = JCAMPDXir.read_jdx_file(raw".\test\tests data\jcamp_python\example_multiline_datasets.jdx")
plot(data.x,data.y)
data.headers

T_data = JCAMPDXir.read_jdx_file(raw".\test\tests data\JCAMP_XYXY.jdx") 
plot(T_data.x,T_data.y)
T_data.headers["NPOINTS"]
length(T_data.y)
head = JCAMPDXir.parse_headers(raw".\test\tests data\JCAMP_XYXY.jdx")
head
include("JCAMPDXir.jl")
jdx_blocks = JCAMPDXir.count_blocks(raw".\test\tests data\jcamp_python\example_compound_file.jdx")
out = JCAMPDXir.read!.(jdx_blocks)
using Plots
plot(out[1].x,out[1].y)
plot!(out[2].x,out[2].y)