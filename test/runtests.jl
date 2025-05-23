using JCAMPDXir,DelimitedFiles,Interpolations
using Test
test_data_folder = joinpath(@__DIR__(),"tests data")
test_file = joinpath(test_data_folder,"JCAMP_test_file.jdx") # this file containes testing spectra in JCAMP-DX specification
no_headers_file = joinpath(test_data_folder,"test_file_no_headers.txt") # data without headers
written_file_name = joinpath(test_data_folder,"written.jdx")
two_column_ascii_file_name = joinpath(test_data_folder,"transmittance.txt")

data_norm(x,y) = sqrt(sum(x->x^2, x .-y))/length(x)
data_norm_rel(x,y) =begin
    dn = data_norm(x,y)
    m = sum(sum(i) for i in zip(x,y))/(2*length(x))
    return dn/m
end
@testset "JCAMPDXir.jl" begin
    # testing by writing and reading the same data
    data = JCAMPDXir.read_jdx_file(test_file) # reading test file
    x_units = data.headers["XUNITS"] # x units of loaded data
    y_units = data.headers["YUNITS"] # y units of loaded data
    JCAMPDXir.write_jdx_file(written_file_name,data.x,data.y,x_units,y_units; yunits =y_units,xunits = x_units ) # writing file with the same units
    data_written = JCAMPDXir.read_jdx_file(written_file_name)
    @show data_norm_rel(data.x,data_written.x)
    @show data_norm_rel(data.y,data_written.y)
    @test data_norm_rel(data.x,data_written.x) <= 1e-5  
    @test data_norm_rel(data.y,data_written.y) <= 1e-5 
    # testing the equality of intermediate data of Int type
    no_headers_data= readdlm(no_headers_file) # reading benchmark data
    x_test = no_headers_data[:,1]
    y_test = Int.(vec(transpose(no_headers_data[:,2:end]))) 
    (x_copy,_,y_int,_) = JCAMPDXir.prepare_jdx_data(data.x,data.y,x_units,y_units; yunits =y_units,xunits = x_units )
    @test data_norm_rel(x_copy[1:8:end],x_test) <= 1e-5  
    @test data_norm_rel(y_int,y_test) <= 1e-5 
    # testing the correctness of units convertion
    T_data = JCAMPDXir.read_jdx_file(two_column_ascii_file_name) # this file stores two-column data without headers
    n_points = 8*div(length(T_data.x),8)
    T_data_x = T_data.x[1:n_points]
    T_data_y = T_data.y[1:n_points]
    @show T_data_x_units = T_data.headers["XUNITS"]
    @show T_data_y_units = T_data.headers["YUNITS"]
    # in T_data x units are MKM and y units are TRANSMISSIVITY
    X_UNITS = values(JCAMPDXir.xNUM2STR) # all possible x units
    Y_UNITS = values(JCAMPDXir.yNUM2STR) # all possible y units
    # now running through all possible units, writing initial T_data to file with converted units
    # reading this file to new data_current and loading values
    # Further we are reading the converted data from file, and checking the values by
    # converting data manually 
    for YU in Y_UNITS
        for XU in X_UNITS
            @show XU
            @show YU
            JCAMPDXir.write_jdx_file(written_file_name,T_data_x,T_data_y,
                                            T_data_x_units,T_data_y_units; 
                                            yunits =YU,xunits = XU,
                                            title="x-units -$(T_data_x_units);yunits-$(T_data_y_units)" ) # writing file with converting units
            data_current = JCAMPDXir.read_jdx_file(written_file_name)
            if XU=="MICROMETERS"
                x_test  = copy(T_data_x)
            elseif XU=="1/CM"
                x_test = 1e4 ./(reverse(T_data_x))
            elseif XU=="NANOMETERS"
                x_test  =  1000 .*T_data_x 
            end

            @show data_norm_rel(x_test,data_current.x)
            @test data_norm_rel(x_test,data_current.x)<1e-1

            if YU == "TRANSMITTANCE" || YU =="REFLECTANCE"
                y_test = copy(T_data_y)

                @show data_norm_rel(y_test,data_current.y)
                @test data_norm_rel(y_test,data_current.y)<1e-1

            end
        end
    end

end