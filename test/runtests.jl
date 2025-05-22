using JCAMPDXir,DelimitedFiles
using Test
test_data_folder = joinpath(@__DIR__(),"tests data")
test_file = joinpath(test_data_folder,"JCAMP_test_file.jdx") # this file containes testing spectra in JCAMP-DX specification
no_headers_file = joinpath(test_data_folder,"test_file_no_headers.txt") # data without headers
written_file_name = joinpath(test_data_folder,"written.jdx")

data_norm(x,y) = sqrt(sum(x->x^2, x .-y))/length(x)
@testset "JCAMPDXir.jl" begin
    # testing by writing and reading the same data
    data = JCAMPDXir.read_jdx_file(test_file) # reading test file
    x_units = data.headers["XUNITS"] # x units of loaded data
    y_units = data.headers["YUNITS"] # y units of loaded data
    JCAMPDXir.write_jdx_file(written_file_name,data.x,data.y,x_units,y_units; yunits =y_units,xunits = x_units ) # writing file with the same units
    data_written = JCAMPDXir.read_jdx_file(written_file_name)
    @test data_norm(data.x,data_written.x) <= 1e-4  
    @test data_norm(data.y,data_written.y) <= 1e-4 
    # testing the equality of intermediate data of Int type
    no_headers_data= readdlm(no_headers_file) # reading benchmark data
    x_test = no_headers_data[:,1]
    y_test = Int.(vec(transpose(no_headers_data[:,2:end]))) 
    (x_copy,_,y_int,_) = JCAMPDXir.prepare_jdx_data(data.x,data.y,x_units,y_units; yunits =y_units,xunits = x_units )
    @test data_norm(x_copy[1:8:end],x_test) <= 1e-3  
    @test data_norm(y_int,y_test) <= 1e-4 
    # testing the correctness of units convertion
    X_UNITS = values(JCAMPDXir.xNUM2STR) # all possible x units
    Y_UNITS = values(JCAMPDXir.yNUM2STR) # all possible y units
    for YU in Y_UNITS
        for XU in X_UNITS
            @show XU
            @show YU
            JCAMPDXir.write_jdx_file(written_file_name,data.x,data.y,x_units,y_units; yunits =YU,xunits = XU )
            data_current = JCAMPDXir.read_jdx_file(written_file_name)
        end
    end

end


data = JCAMPDXir.read_jdx_file(test_file) # reading test file
x_units = data.headers["XUNITS"] # x units of loaded data
y_units = data.headers["YUNITS"] # y units of loaded data
JCAMPDXir.write_jdx_file(written_file_name,data.x,data.y,x_units,y_units; yunits ="ARBITRARY UNITS",xunits = "NANOMETERS" )
x_copy = copy(data.x)
xconvert!(x_copy[2:end],"1/CM","NANOMETERS")
x_copy