const xSTR2NUM = Dict("MKM"=>1,"MICROMETERS"=>1,"WAVELENGTH (UM))"=>1, #('micrometers', 'um', 'wavelength (um)')
"1/CM"=>2,"CM-1"=>2,"CM^-1"=>2, # ('1/cm', 'cm-1', 'cm^-1')
"NM"=>3,"NANOMETERS"=>3,"WAVELENGTH (NM))"=>3)     
const xNUM2STR = Dict(1=>"MICROMETERS",2=>"1/CM",3=>"NANOMETERS")
const ySTR2NUM = Dict("TRANSMITTANCE"=>1,"T"=>1,"REFLECTANCE"=>2,"R"=>2,
        "ABSORBANCE"=>3,"A"=>3,"KUBELKA-MUNK"=>4,"ARBITRARY UNITS"=>5,"A.U."=>5)
const yNUM2STR = Dict(1=>"TRANSMITTANCE",2=>"REFLECTANCE",3=>"ABSORBANCE",4=>"KUBELKA-MUNK",5=>"ARBITRARY UNITS") 
const supported_x_units = join(keys(xSTR2NUM),",");
const supported_y_units = join(keys(ySTR2NUM),",");
#const MAX_POINTS_IN_LINE = 15 # this limits the number of point per single line 
const SUPPORTED_JCAMPDX_VERSION =4.24
abstract type sUnits{T} end
struct yUnits{T}<:sUnits{T}
yUnits(u::String) = begin
return haskey(ySTR2NUM,uppercase(u)) ? new{ySTR2NUM[uppercase(u)]}() : new{5}()
end
end
struct xUnits{T} <:sUnits{T}
xUnits(u::String) = begin
return haskey(xSTR2NUM,uppercase(u)) ? new{xSTR2NUM[uppercase(u)]}() : error("this x units are not supported, possible units are: "*supported_x_units)
end
end
units(::xUnits{T}) where T = xNUM2STR[T]
units(::yUnits{T}) where T = yNUM2STR[T]
"""
    convert!(x::AbstractArray,::sUnits{T},::sUnits{T}) where T

This functions can be used convert y and x units 

supported x-units names: $(supported_x_units)

supported y-units names: $(supported_y_units)

# Example
```
julia> convert!([1,2,3],xUnits("MKM"),xUnits("1/cm"))) 
```
"""
    convert!(x::AbstractArray,::sUnits{T},::sUnits{T}) where T = x 
    convert!(x::AbstractArray,::sUnits{T},::sUnits{P}) where T where P = x 
#x - converters
    convert!(x::AbstractArray,::xUnits{1},::xUnits{2}) = @. x = 1e4/x #mkm=>1/cm
    convert!(x::AbstractArray,::xUnits{2},::xUnits{1}) = @. x = 1e4/x #1/cm=>mkm
    convert!(x::AbstractArray,::xUnits{1},::xUnits{3}) = @. x = 1e3*x #mkm=>nm
    convert!(x::AbstractArray,::xUnits{3},::xUnits{1}) = @. x = 1e-3*x #nm=>mkm
    convert!(x::AbstractArray,::xUnits{2},::xUnits{3}) = @. x = 1e7/x #1/cm=>nm
    convert!(x::AbstractArray,::xUnits{3},::xUnits{2}) = @. x = 10.0/x #nm=>1/cm
    
# Y-converters
    convert!(x::AbstractArray,::yUnits{1},::yUnits{3}) = @. x=-log10(x) #T->A
    convert!(x::AbstractArray,::yUnits{3},::yUnits{1}) = @. x=10^(-x)# A->T
    convert!(x::AbstractArray,::yUnits{2},::yUnits{3}) = @. x=-log10(x)#R->A
    convert!(x::AbstractArray,::yUnits{3},::yUnits{2}) = @. x=10^(-x) #A->R   
    convert!(x::AbstractArray,::yUnits{1},::yUnits{4}) = @. x= (1 - x^2)/(2*x) #T->K-M
    convert!(x::AbstractArray,::yUnits{4},::yUnits{1}) = @. x=-x  + sqrt(x^2 + 1)# K-M->T
    convert!(x::AbstractArray,::yUnits{2},::yUnits{4}) = @. x=(1 - x^2)/(2*x)#R->K-M
    convert!(x::AbstractArray,::yUnits{4},::yUnits{2}) = @. x= -x  + sqrt(x^2 + 1) #K-M->R    
    
    """
    xconvert!(x::AbstractArray,input_units::String,output_units::String)

    Converts the values of `x` from `input_units` to `output_units`.
    supported x-units names: $(supported_x_units)

# Example
```
julia> xconvert!([1,2,3],"MKM",xUnits"1/cm")) 
```
"""
xconvert!(x::AbstractArray,input_units::String,output_units::String) = convert!(x,xUnits(input_units),xUnits(output_units))
    """
    yconvert!(y::AbstractArray,input_units::String,output_units::String)

Converts the values of `y` from `input_units` to `output_units`.
supported x-units names: $(supported_x_units)
All units can be written both in lower- and in uppercase, `T`,`R` and `A` stay for 
a shorthand for `TRANSMITTANCE`,`REFLECTANCE` and `ABSORBANCE`

# Example
```julia
julia> yconvert!([1,2,3],"R","KUBELKA-MUNK")) 
```
"""
yconvert!(y::AbstractArray,input_units::String,output_units::String) = convert!(y,yUnits(input_units),yUnits(output_units))
    
    """
    write_jdx_file(file_name,x::Vector{Float64},y::Vector{Float64},x_units::String="1/CM",
        y_units::String="TRANSMITTANCE"; kwargs...)

Saves infrared spectrum given as two vectors `x` and `y` to JCAMP-DX file `file_name`.
Input vector should be of the same size. Currently the package suppots only `(X++(Y..Y))`
table data format. `JCAMP-DX 4.24` demands 80 symbols per file line of Y-data and 88 total symbols per line.
The y-vector is stored in eight columns, thus the total number of points should be a multiple of eight.
All last `mod(length(y),8)` points of y dtaa will not be written to the file. It is preferable that all x-data is
sorted in ascending order and spaced uniformly. If it is not the case, or if the units conversion 
is envolved, function will automatically interpolate and sort the data on uniformly spaced grid.

x_units  - units of x data, must be one of  $(supported_x_units)

y_units -  units of y data, must be one of $(supported_y_units) 

Further any keword arguments can be provided, all of them will be written to the head of the file.
All keyword arguments appear in the file in uppercase.

Most impostant keywords are 

    TITLE - the title of the file (it is always on top of the file)
    XUNITS - x data units saved to file,  must be one of  $(supported_x_units)
    YUNITS - y data units saved to file, must be one of $(supported_y_units) 

If `x_units` (function's fourth argument) are not equal to the key-word argument XUNITS than the function converts x-values before saving to file see [`xconvert!`](@ref)

If `y_units` (function's fifth argument) are not equal to the key-word argument XUNITS than the function converts y-values before saving to file see [`yconvert!`](@ref)

# Example
```julia
julia> using JCAMPDXir
julia> filename = joinpath(@__DIR__,"test.jdx")
julia> write_jdx_file(filename,[1,2,3,4,5,6,7,8],rand(8),"MKM","T",title = "new file",XUNIT="1/CM",YUNITS="KUBELKA-MUNK") 

```
"""
function write_jdx_file(file_name,x::Vector{Float64},y::Vector{Float64},x_units::String="1/CM",
        y_units::String="TRANSMITTANCE"; kwargs...)

        (x_copy,_,y_int,headers) = prepare_jdx_data(x,y,x_units,y_units; kwargs...)
        y_columns_number = round(Int,80/(ndigits(JCAMPDXir.YMAX_INT)+1)) # number of columns of y-data
        npoints = length(y_int) # total number of data points
        headers["NPOINTS"] = Float64(y_columns_number*round(Int,npoints/y_columns_number))
        fmt = Printf.Format(" %d"^y_columns_number*"\r\n")
        open(file_name,"w", lock = true) do io
            for (k,v) in headers
                k == "NPOINTS" ? v_str = string(round(Int,v)) : v_str=string(v)
                line = "##"*k*"="*v_str*"\r\n"
                write(io,line)
            end
            counter = 1
            line_index = 1
            while counter <= npoints
                start_ind = counter
                end_ind = counter+y_columns_number-1
                if end_ind<=npoints
                    line = Printf.format(fmt,y_int[start_ind:end_ind]...)
                    x_str = @sprintf("%.8f",x_copy[start_ind])
                    remained_length = 89 - length(line)
                    remained_length>length(x_str) ? write(io,x_str*line) : write(io,x_str[1:remained_length]*line)
                else
                    break
                end
                counter = end_ind+1
                line_index+=1
            end
            write(io,"##END=\r\n")
        end
        return nothing
    end
    """
    prepare_jdx_data(x::Vector{Float64},y::Vector{Float64},x_units::String="1/CM",
                                                    y_units::String="TRANSMITTANCE"; kwargs...)

Prepares data to be written using [`write_jdx_file`](@ref)
"""
function prepare_jdx_data(x::Vector{Float64},y::Vector{Float64},x_units::String="1/CM",
                                                    y_units::String="TRANSMITTANCE"; kwargs...)
        @assert length(x)==length(y)
        x_init = copy(x)
        y_init = copy(y)
        headers = copy(default_headers)
        for (k,v) in kwargs
            k_str = uppercase(string(k))
            isa(v,String) ? v=uppercase(v) : nothing
            headers[k_str] = v 
        end
        #must ensure the the  "XYDATA" element of headers is the last because they are writen in the order
        v = pop!(headers,"XYDATA")
        push!(headers,"XYDATA"=>v)
        xconvert!(x_init, x_units, headers["XUNITS"])
        yconvert!(y_init, y_units, headers["YUNITS"])
        # after unit conversion there can be the case that some data is
        # NaN or Inf, thus need the isfinite check
        inf_inds = Vector{Int}([])
        for i in eachindex(x_init)
            !isfinite(x_init[i]) || !isfinite(y_init[i]) ? push!(inf_inds,i) : nothing
        end 
        if !isempty(inf_inds)
            deleteat!(x_init,inf_inds)
            deleteat!(y_init,inf_inds)
        end
        x_factor = headers["XFACTOR"] 
        n_points = length(x_init)
        if is_linspaced(x_init)# checks if all coordinates are equally spaced, if not - performing interpolation
            x_copy   = x_init
            y_copy = y_init
            x_copy ./= x_factor
            if !issorted(x_copy)
                y_int = sortperm(x_copy)
                @. x_copy=x_copy[y_int]
                @. y_copy=y_copy[y_int]
            else
                y_int = Vector{Int}(undef,n_points) 
            end
        else # if x is not equally spaced we perform linear interpolation
            x_copy = collect(range(minimum(x_init),maximum(x_init),n_points)) 
            if !issorted(x_init)
                y_int = sortperm(x_init)
                interpolator =  linear_interpolation(x_init[y_int],y_init[y_int])
            else
                y_int = Vector{Int}(undef,n_points)
                interpolator =  linear_interpolation(x_init,y_init)
            end
            y_copy = interpolator(x_copy)
            x_copy./=x_factor
            interpolator = nothing
        end
        cur_date_time = string(now())
        ind = findfirst("T",cur_date_time)[1]
        headers["NPOINTS"] = Float64(n_points)
        headers["DATE"] = cur_date_time[1:ind-1]
        headers["TIME"] = cur_date_time[ind+1:end]
        headers["FIRSTX"] = x_copy[begin]
        headers["FIRSTY"] = y_copy[begin]
        (headers["MINY"],headers["MAXY"]) = extrema(y_copy)
        (headers["MINX"],headers["MAXX"]) = extrema(x_copy)
        headers["FIRSTX"] = headers["MINX"]
        headers["LASTX"] = headers["MAXX"]
        y_factor = (headers["MAXY"]/YMAX_INT)
        @. y_int = round(Int,y_copy/y_factor) #filling integer values
        headers["YFACTOR"] = y_factor
        return (x_copy,y_copy,y_int,headers)
    end
    function is_linspaced(x::Vector{T}) where T<:Number
        if length(x)<=2
            return true
        end
        dx1 = x[2] - x[1]
        for i in eachindex(x)[2:end-1]
           dx = x[i + 1] - x[i]
           isapprox(dx,dx1,rtol=1e-8) ? continue : return false
        end
        return true
    end