include("JCAMPDXir.jl")
test_data_folder = joinpath(@__DIR__(),"..")
test_data_folder = joinpath(test_data_folder,"test","tests data")
test_file = joinpath(test_data_folder,"JCAMP_test_file.jdx") 
using .JCAMPDXir
jdx = JDXfile(test_file)
@code_warntype JCAMPDXir.read!(jdx) # 
data = JCAMPDXir.read!(jdx) # reading test file
x_units = data.headers["XUNITS"] # x units of loaded data
y_units = data.headers["YUNITS"] # y units of loaded data  


# 
using StaticArrays,BenchmarkTools
function add1(y_data, current_line::String,number_of_y_point_per_chunk::Int)
    #SVector{Float64,}
   #data_intermediate = map((X)->Base.parse(Float64,X), split(current_line))
   data_intermediate = SVector{number_of_y_point_per_chunk+1}(Base.parse.(Float64,split(current_line)))
   starting_index =1
   y_data[starting_index:starting_index+number_of_y_point_per_chunk-1] = data_intermediate[2:end]
   return data_intermediate[1]::Float64 # returns x-value for checks
end

function add0(y_data, current_line::String,number_of_y_point_per_chunk::Int)
    #SVector{Float64,}
   #data_intermediate = map((X)->Base.parse(Float64,X), split(current_line))
   data_intermediate = Base.parse.(Float64,split(current_line))
   starting_index =1
   y_data[starting_index:starting_index+number_of_y_point_per_chunk-1] = data_intermediate[2:end]
   return data_intermediate[1]::Float64 # returns x-value for checks
end
function add2(y_data, current_line::String,number_of_y_point_per_chunk::Int)
    #SVector{Float64,}
   #data_intermediate = map((X)->Base.parse(Float64,X), split(current_line))
   data_intermediate = Vector{Float64}(undef,number_of_y_point_per_chunk+1)
   for (i,s) in enumerate(eachsplit(current_line))
        data_intermediate[i]= Base.parse(Float64,s)
   end
   #data_intermediate = Base.parse.(Float64,split(current_line))
   starting_index =1
   y_data[starting_index:starting_index+number_of_y_point_per_chunk-1] = data_intermediate[2:end]
   return data_intermediate[1]::Float64 # returns x-value for checks
end
function add3(y_data, current_line::String,number_of_y_point_per_chunk::Int)
    #SVector{Float64,}
   #data_intermediate = map((X)->Base.parse(Float64,X), split(current_line))
   data_intermediate = MVector{number_of_y_point_per_chunk+1,Float64}(undef)
   for (i,s) in enumerate(eachsplit(current_line))
        data_intermediate[i]= Base.parse(Float64,s)
   end
   #data_intermediate = Base.parse.(Float64,split(current_line))
   starting_index =1
   y_data[starting_index:starting_index+number_of_y_point_per_chunk-1] = data_intermediate[2:end]
   return data_intermediate[1]::Float64 # returns x-value for checks
end

function add4(y_data, current_line::String,number_of_y_point_per_chunk::Int)
    #SVector{Float64,}
   #data_intermediate = map((X)->Base.parse(Float64,X), split(current_line))
   data_intermediate = Vector{Float64}(undef,number_of_y_point_per_chunk+1)
   for (i,s) in enumerate(eachsplit(current_line))
         v = Base.tryparse(Float64,s)
         isnothing(v) ? data_intermediate[i]=-1 : data_intermediate[i]= v
   end
   #data_intermediate = Base.parse.(Float64,split(current_line))
   starting_index =1
   y_data[starting_index:starting_index+number_of_y_point_per_chunk-1] = data_intermediate[2:end]
   return data_intermediate[1]::Float64 # returns x-value for checks
end
y_data = Vector{Float64}(undef,8)
@benchmark add1(y_data,"15794.7 122647 117553 111281 147129 161624 182712 353010 764665",8)
@benchmark add0(y_data,"15794.7 122647 117553 111281 147129 161624 182712 353010 764665",8)
@benchmark add2(y_data,"15794.7 122647 117553 111281 147129 161624 182712 353010 764665",8)
@benchmark add3(y_data,"15794.7 122647 117553 111281 147129 161624 182712 353010 764665",8)
@benchmark add4(y_data,"15794.7 122647 117553 111281 147129 161624 182712 353010 764665",8)


######################
include("JCAMPDXir.jl")
test_data_folder = joinpath(@__DIR__(),"..")
test_data_folder = joinpath(test_data_folder,"test","tests data")
test_file = joinpath(test_data_folder,"JCAMP_test_file.jdx") # this file containes testing spectra in JCAMP-DX specification
no_headers_file = joinpath(test_data_folder,"test_file_no_headers.txt") # data without headers
written_file_name = joinpath(test_data_folder,"written.jdx")
two_column_ascii_file_name = joinpath(test_data_folder,"transmittance.txt")
T_data = JCAMPDXir.read_jdx_file(two_column_ascii_file_name,delimiter="\t")