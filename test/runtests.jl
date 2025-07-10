
using Test
#using JCAMPDXir
using GitHub,Downloads,DelimitedFiles,Interpolations

test_data_folder = joinpath(@__DIR__(),"tests data")
python_package_test_data = joinpath(test_data_folder,"jcamp_python")
files_from_jcamp_py = filter(r->occursin(".jdx",r), readdir(python_package_test_data))
#cheking files from jcamp repository for python
try
    owner = "nzhagen"
    repo = "jcamp"
    path = "data/infrared_spectra"
    print("Refreshing testing data from ($(joinpath(owner,repo)))...")
    R = GitHub.repo(joinpath(owner,repo))
    files, _ = GitHub.directory(R, path)
    
    for f in files
        fname = f.name
        !occursin(".jdx",fname) || fname ∈ files_from_jcamp_py ? continue : nothing
        Downloads.download(string(f.download_url),joinpath(python_package_test_data,fname))
    end
    println("ok")
catch 
    println("unable to read data from repo")
end

test_file = joinpath(test_data_folder,"JCAMP_test_file.jdx") # this file containes testing spectra in JCAMP-DX specification
no_headers_file = joinpath(test_data_folder,"test_file_no_headers.txt") # data without headers
written_file_name = joinpath(test_data_folder,"written.jdx")
two_column_ascii_file_name = joinpath(test_data_folder,"transmittance.txt")
XYXY_test_file = joinpath(test_data_folder,"JCAMP_XYXY.jdx")
JCAMPDX5_test_file = joinpath(test_data_folder,"JCAMPDX5.jdx")

JCAMPDX_ASDF_file = joinpath(python_package_test_data,"isopropanol_ASDF.jdx")
ASDF_file_decoded = joinpath(test_data_folder,"isopropanol_ASDF_decoded.txt") # this file was decoded using jcamp package 


data_norm(x,y) = sqrt(sum(x->x^2, x .-y))/length(x)
data_norm_rel(x,y) =begin
    dn = data_norm(x,y)
    m = sum(sum(i) for i in zip(x,y))/(2*length(x))
    return dn/m
end
A2T(x) = begin # function to convert the absorbance to transmission 
    return  @. 10^(-x) #A = log10(I0/I) = -log10(T) => T=10^(-A)
end
KM2T(x) = begin # function to convert Kubelka-Munk to transmittance
    return  @. 10^(-x) #A = log10(I0/I) = -log10(T) => T=10^(-A)
end

@testset "JCAMPDXir.jl" begin
    println("____________________")
    println("\ntesting files from $(python_package_test_data)" )
    for f in files_from_jcamp_py
       try 
            print("reading $(f) ")
            JCAMPDXir.read_jdx_file(joinpath(python_package_test_data,f))
            @test true
            print("...ok \n")
       catch ex
            print("...reading error \n")
            @show ex
            @test false
       end
    end
    println("____________________")
    println("\nReading XYYdata")
    data = JCAMPDXir.read_jdx_file(test_file) 
    for (i,y) in enumerate(0.0523676*[1811180; 1261721; 2334144; 2496642;
                                    2864863; 3238081; 3165312; 3090651; 
                                    2718186; 2951604; 3192506; 3573356;
                                    3591102; 3219943; 3321428; 3797139])
        @test data.y[i] ≈ y
    end
    println("____________________")    
    println("\nReading XYXYdata")
    data = JCAMPDXir.read_jdx_file(XYXY_test_file) 
    for (i,y) in enumerate(eachrow([11995.21 32112;
                              11991.36 32505;
                              11987.5 32727;
                              11983.64 33481;
                              11979.79 33798;
                              11975.93 33076;
                              11972.07 32550;
                              11968.22 33189]))
       @test data.y[i] ≈ 1e-6*y[2]
       @test data.x[i] ≈ y[1]   
    end 
    println("____________________")    
    println("\nReading JCAMP-DX-5.01")
    data = JCAMPDXir.read_jdx_file(JCAMPDX5_test_file)
    for (i,y) in enumerate(vec([114292232 112522752 111130888 110025904 109160424 108469448 107820264 107003944 105776976 103944576 101429344 98357408 94981288]))
        @test data.y[i] ≈ 1.000000E-009*y 
    end


    println("____________________")    
    println("\nReading ASDF file format")
    data = JCAMPDXir.read_jdx_file(JCAMPDX_ASDF_file)
    asdf_data = readdlm(ASDF_file_decoded) # reading data decoded using python jcamp parser
    @test begin
        discr = 0
        for (y1,y2) in zip(asdf_data[:,2],data.y)
            discr +=abs(y1-y2)
        end
        discr/length(data.y)
    end ≈ 0

    println("____________________")
    println("\nTesting by writing and reading the same data")
    data = JCAMPDXir.read_jdx_file(test_file) # reading test file
    x_units = data.headers["XUNITS"] # x units of loaded data
    y_units = data.headers["YUNITS"] # y units of loaded data
    JCAMPDXir.write_jdx_file(written_file_name,data.x,data.y,x_units,y_units; yunits =y_units,xunits = x_units ) # writing file with the same units
    data_written = JCAMPDXir.read_jdx_file(written_file_name)
    dx_norm = data_norm_rel(data.x,data_written.x)
    dy_norm = data_norm_rel(data.y,data_written.y)
    println("relative error: S(x) =  $(dx_norm), S(y) = $(dy_norm) ")
    @test dx_norm<= 1e-5  
    @test dy_norm<= 1e-5 
    println("____________________")
    println("\nTesting the equality of intermediate data of integer type used to compress the y-data")
    no_headers_data= readdlm(no_headers_file) # reading benchmark data
    x_test = no_headers_data[:,1]
    y_test = Int.(vec(transpose(no_headers_data[:,2:end]))) 
    (x_copy,_,y_int,_) = JCAMPDXir.prepare_jdx_data(data.x,data.y,x_units,y_units; yunits =y_units,xunits = x_units )
    dx_norm = data_norm_rel(x_copy[1:8:end],x_test)
    dy_norm = data_norm_rel(y_int,y_test)
    println("relative error: S(x) =  $(dx_norm), S(y) = $(dy_norm) ")
    @test dx_norm <= 1e-5  
    @test dy_norm <= 1e-5 
    println("____________________")
    println("\nTesting the correctness of units convertion ")
    T_data = JCAMPDXir.read_jdx_file(two_column_ascii_file_name,delimiter="\t") 
    # this file stores two-column data with XUNITS and YUNITS the data is delimited with tab
    n_points = 8*div(length(T_data.x),8) # the number of points should be an integer multiplyer of 8
    T_data_x = T_data.x[1:n_points]
    T_data_y = T_data.y[1:n_points]
    T_data_x_units = T_data.headers["XUNITS"]
    T_data_y_units = T_data.headers["YUNITS"]
    # in T_data x units are MKM and y units are TRANSMISSIVITY
    X_UNITS = values(JCAMPDXir.xNUM2STR) # all possible x units
    Y_UNITS = values(JCAMPDXir.yNUM2STR) # all possible y units
    # now running through all possible units, writing initial T_data to file with converted units
    # reading this file to new data_current and loading values
    # Further we are reading the converted data from file, and checking the values by
    # converting data manually 
    for YU in Y_UNITS
        for XU in X_UNITS
            println("x : $(T_data_x_units) => $(XU) , y :  $(T_data_y_units) => $(YU)")
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

            dx_norm = data_norm_rel(x_test,data_current.x)
            @test dx_norm < 1e-1

            if YU == "TRANSMITTANCE" || YU =="REFLECTANCE" || YU=="ARBITRARY UNITS"# units are the same as the reference data
                interp1= linear_interpolation(x_test,T_data_y)
            elseif YU == "ABSORBANCE"
                interp1= linear_interpolation(x_test,A2T(T_data_y))
            elseif YU == "KUBELKA-MUNK"
                interp1= linear_interpolation(x_test,KM2T(T_data_y))
            end
            y_test = interp1(data_current.x)
           # @show length(y_test) length(data_current.y)
           dy_norm = data_norm_rel(y_test,data_current.y)
           println("relative error: S(x) =  $(dx_norm), S(y) = $(dy_norm) \n")
            @test dy_norm < 1e-2
        end
    end

end