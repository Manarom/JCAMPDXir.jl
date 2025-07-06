#debug script
using Revise, Plots,OrderedCollections,StaticArrays,DelimitedFiles
test_data_folder = joinpath(".","test","tests data")
python_package_test_data = joinpath(test_data_folder,"jcamp_python")
ASDF_file_decoded = joinpath(test_data_folder,"isopropanol_ASDF_decoded.txt")
asd_file = joinpath(python_package_test_data,"isopropanol_ASDF.jdx")
asdf_data = readdlm(ASDF_file_decoded)
includet("JCAMPDXir.jl")
plot(asdf_data[:,1],asdf_data[:,2])
asdf_line = "2561A058O445j574k394k500J085J232K650J052o294P237k468p57m255K018L444j210"
asdf_data[1]
asdf_headers = JCAMPDXir.parse_headers(asd_file)
y_asdf = vec(asdf_data[:,2])/asdf_headers["YFACTOR"]
x_asdf = vec(asdf_data[:,1])/asdf_headers["XFACTOR"]
x_asdf[1]
dif_keys = keys(JCAMPDXir.DIF_digits)
ln = replace(asdf_line,[k=>" "*string(k) for k in keys(JCAMPDXir.DIF_digits)]...)
y_asdf[1] 
y_asdf[2]
y_asdf[3]
ln = asdf_line |> JCAMPDXir.SQZ |> JCAMPDXir.DIF
ln = replace(ln,JCAMPDXir.DIF_digits...)
db = [parse(Float64,s) for s in eachsplit(ln)]
dbv = @view db[2:end] 
dbv3 = @view db[3:end]
for i in eachindex(db)[3:end]
    db[i] = db[i-1] + db[i]
end
db
using LinearAlgebra
norm(db[2:end] .- y_asdf[1:17])
y_asdf[3]
ln = JCAMPDXir.SQZ(ln)
replace(ln,)


d = JCAMPDXir.DataBuffer(MVector{3,Float64}(undef))
JCAMPDXir.fill_data_buffer!(d,"2.3 4.6 7.8",isspace)
d.buffer


d2= JCAMPDXir.DataBuffer()
line = "575.17-3042244 1597332-970474 1254921-1092859 2206136"
jdx = JCAMPDXir.JDXblock()
JCAMPDXir.addline!(jdx,d2,line)
jdx.y_data
data = JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,"benzene.jdx")) # reading test file
plot(data.x,data.y)
a = JCAMPDXir.DataBuffer()
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

T_data = JCAMPDXir.read_jdx_file(raw".\test\tests data\JCAMPDX5.jdx") 
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