#debug script
using Revise
using Plots,OrderedCollections
test_data_folder = joinpath(".","test","tests data")
python_package_test_data = joinpath(test_data_folder,"jcamp_python")
#includet("JCAMPDXir.jl")
include("JCAMPDXir.jl")
using .JCAMPDXir
line = "575.17-3042244 1597332-970474 1254921-1092859 2206136"
a= Float64[]
JCAMPDXir.fill_data_chunk!(a,"575.17 3042244 1597332-970474 1254921  1092859 2206136",isspace)
a


using StaticArrays
v = MVector{7,Float64}(a)
JCAMPDXir.fill_data_chunk!(v,line,isspace)

v
jdx = JCAMPDXir.JDXblock()
JCAMPDXir.addline!(JCAMPDXir.XYYline,jdx,
                                line)
jdx.y_data
begin
    data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,"1-butanol.jdx")) # reading test file
    return plot(data.x,data.y)
end
a = [1.0; 1.0; 1.0; 1.0; 1.0]
JCAMPDXir.split_PAC_string!(a,3,"+345-23.0-4")
a

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
        error_counter +=1
        push!(problem_files_index,error_counter=>Pair(i,msg))
    end
end
problem_files_index[1]
error_counter
println("error $(100*error_counter/length(files_from_jcamp_py)) %")
files_from_jcamp_py[2]
data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,files_from_jcamp_py[2])) # reading test file
plot(data.x,data.y)
data.y[end]

data_headers = JCAMPDXir.parse_headers(joinpath(python_package_test_data,files_from_jcamp_py[36]))

data = JCAMPDXir.read_jdx_file(raw".\test\tests data\JCAMP_test_file.jdx")
plot(data.x,data.y)
data.headers

T_data = JCAMPDXir.read_jdx_file(raw".\test\tests data\jcamp_python\example_compound_file.jdx") 
plot(T_data.x,T_data.y)
data.headers["NPOINTS"]
