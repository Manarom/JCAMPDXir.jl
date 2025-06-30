
# module to read JCAMP-DX=4.24 file formats

module JCAMPDXir
using Dates,Interpolations,OrderedCollections,Printf,StaticArrays
export JDXblock,
            read!,
            read_jdx_file,
            write_jdx_file,
            parse_headers,
            xconvert!,
            yconvert!,
            addline!,
            XYYline
"""
This JCAMP-DX (infrared) file format was developed for the exchange of infrared spectra between different laboratories.
For general description of format refer to  [`UIPAC.ORG (pdf file)`](https://iupac.org/wp-content/uploads/2021/08/JCAMP-DX_IR_1988.pdf)
In addition to the spectra themselves, the file also stores metadata containing information about the units of measurement
and the conditions under which the spectra were acquired. A detailed specification of the format is provided via the link. 
This package was written for a specific task and cannot be considered as a complete implementation of all JCAMP-DX format specification; 
however, it implements some of the basic functionality to read [`read_jdx_file`](@ref) and write [`write_jdx_file`](@ref) JCAMP-DX files
according to the `JCAMP-DX=4.24`.

JCAMP file content example: 

            ##TITLE=1 
            ##JCAMP-DX=4.24
            ##DATATYPE=INFRARED SPECTRUM
            ##DATE=2021/10/17
            ##TIME=11:30
            ##XUNITS=1/CM
            ##YUNITS=TRANSMITTANCE
            ##XFACTOR=1
            ##YFACTOR=0.00699183
            ##FIRSTX=0
            ##LASTX=15801.4671743
            ##FIRSTY=11746.3412893072
            ##MAXX=15801.4671743
            ##MINX=0
            ##MAXY=3.75371e+006
            ##MINY=4040.25
            ##NPOINTS=16384
            ##XYDATA=(X++(Y..Y))
            0 1680010 821286 2148133 1505245 1537124 1367661 1147725 1134981
            7.71603 1166853 1213186 1029828 1067595 1135904 1195128 1157134 1150556
            15.4321 1266743 1164401 1014224 1022338 999780 1138781 1208950 1161258
            .
            .
            ##END=


"""
JCAMPDXir
    const YMAX_INT = round(Int64,typemax(Int32)/4)-1
    const default_headers = OrderedDict{String,Union{Float64,String}}(
            "TITLE"=>"NO TITLE",
            "JCAMP-DX"=>4.24,
            "DATATYPE"=>"INFRARED SPECTRUM",
            "DATE"=>"2021/10/17",
            "TIME"=>"11:30",
            "XUNITS"=>"1/CM",
            "YUNITS"=>"TRANSMITTANCE",
            "XFACTOR"=>1.0,
            "YFACTOR"=>0.00699183,
            "FIRSTX"=>0.0,
            "LASTX"=>15801.4671743,
            "FIRSTY"=>11746.3412893072,
            "MAXX"=>15801.4671743,
            "MINX"=>0.0,
            "MAXY"=>3.75371e+006,
            "MINY"=>4040.25,
            "NPOINTS"=>16384.0,
            "XYDATA"=>"(X++(Y..Y))"
    )
    const xSTR2NUM = Dict("MKM"=>1,"MICROMETERS"=>1,"WAVELENGTH (UM))"=>1, #('micrometers', 'um', 'wavelength (um)')
                        "1/CM"=>2,"CM-1"=>2,"CM^-1"=>2, # ('1/cm', 'cm-1', 'cm^-1')
                        "NM"=>3,"NANOMETERS"=>3,"WAVELENGTH (NM))"=>3)     
    const xNUM2STR = Dict(1=>"MICROMETERS",2=>"1/CM",3=>"NANOMETERS")
    const ySTR2NUM = Dict("TRANSMITTANCE"=>1,"T"=>1,"REFLECTANCE"=>2,"R"=>2,
                                "ABSORBANCE"=>3,"A"=>3,"KUBELKA-MUNK"=>4,"ARBITRARY UNITS"=>5,"A.U."=>5)
    const yNUM2STR = Dict(1=>"TRANSMITTANCE",2=>"REFLECTANCE",3=>"ABSORBANCE",4=>"KUBELKA-MUNK",5=>"ARBITRARY UNITS") 
    const supported_x_units = join(keys(xSTR2NUM),",");
    const supported_y_units = join(keys(ySTR2NUM),",");
    const X_VIOLATION_CRITERIA = 1e-2
    const DELTAX_VIOLATION_CRITERIUM = 1e-3
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
    # strange symbols
    """
    Squezzed from digits dictionary

    """
    const SQZ_digits = Dict('@'=>"+0", 'A'=>"+1", 'B'=>"+2", 
                            'C'=>"+3", 'D'=>"+4",'E'=>"+5",
                            'F'=>"+6", 'G'=>"+7", 'H'=>"+8", 
                            'I'=>"+9",'a'=>"-1", 'b'=>"-2", 
                            'c'=>"-3", 'd'=>"-4", 'e'=>"-5", 
                            'f'=>"-6", 'g'=>"-7", 'h'=>"-8", 
                            'i'=>"-9")#','=>" " '+'=>"+",  '-'=>"-",  
                            #','=>" ") misleading with XY...XY delimiter and PAC form
	const DIF_digits = Dict('%'=>"+0", 'J'=>"+1",  'K'=>"+2",  
                            'L'=>"+3",  'M'=>"+4",  'N'=>"+5",  
                            'O'=>"+6",  'P'=>"+7",  'Q'=>"+8",
                            'R'=>"+9",  'j'=>"-1", 'k'=>"-2", 
                            'l'=>"-3", 'm'=>"-4", 'n'=>"-5",
                            'o'=>"-6", 'p'=>"-7", 'q'=>"-8",
                            'r'=>"-9")

    const DUP_digits = Dict('S'=>1, 'T'=>2, 'U'=>3, 'V'=>4, 
                            'W'=>5, 'X'=>6, 'Y'=>7, 'Z'=>8,
                             's'=>9)
    
    const SQZ_regexp = Regex("["*join(keys(SQZ_digits))*"]" )
    const DIF_regexp = Regex("["*join(keys(DIF_digits))*"]" )
    const DUP_regexp = Regex("["*join(keys(DUP_digits))*"]" )

    """
Types used to decode the data line before parsing to digits 
By default there is no decoding `NoDecoding`, but if the first line of data 
in block contains `SQZ_digits` than it marks entire block as `SQZDecoding`
by setting the decoding field of `JDXblock` object. String decoding is performed 
right after the reading the line of data from file.

    """
    abstract type Decoding end
    struct SQZDecoding<:Decoding end
    struct NoDecoding<:Decoding end
    struct DIFDecoding<:Decoding end
    struct DUPdecoding<:Decoding end
    struct DIFDUPDecoing<:Decoding end

    decode(::Type{NoDecoding},s::AbstractString) = s
    decode(::Type{SQZDecoding},s) = replace(s,SQZ_digits...)
    decode(::Type{DIFDecoding},s) = replace(s,DIF_digits...)
    (::Type{T})(s::AbstractString) where T<:Decoding = decode(T,s)

    """
    is_PAC_string(s::AbstractString)

Checks if string is of PAC format (+ or - signs are used to separate digits)
"""
is_PAC_string(s::AbstractString) =occursin('-',s)||occursin('+',s)
"""
    is_SQZ_string(s::AbstractString)

Checks if the string is coded using SQZ fomat
"""
is_SQZ_string(s::AbstractString) = occursin(SQZ_regexp,s)
is_DIF_string(s::AbstractString) = occursin(DIF_regexp,s)
is_DUP_string(s::AbstractString) = occursin(DUP_regexp,s)

function get_decoding_type(s::AbstractString)
    is_SQZ_string(s) && return SQZDecoding
    is_DIF_string(s) && is_DUP_string(s) && return DIFDUPDecoing
    is_DIF_string(s) && return DIFDecoding
    is_DUP_string(s) && return DUPdecoding
    return NoDecoding
end

    abstract type DATAline end
    struct XYYline<:DATAline end
    struct XYXYline<:DATAline end

    macro tuple_unpack(N::Int, x)
        @assert N >= 1
        expr = Expr(:tuple)
        push!(expr.args, quote
          begin
            (val, state) = iterate($(esc(x)))
            val
          end
        end)
        for i = 2:N
          push!(expr.args, quote
            begin
              (val, state) = iterate($(esc(x)), state)
              val
            end
          end)
        end
        expr
    end
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
    const IntOrNothing = Union{Int,Nothing}
    """
    Stored parsed data
Must be filled using [`read!`](@ref) function
"""
    mutable struct JDXblock
        # Main struct, prepares  file name, parses data and data headers
        file_name::String
        data_headers::Dict{String,Union{String,Float64}}
        x_data::Vector{Float64}
        y_data::Vector{Float64}
        is_violated_flag::Vector{Bool} # 
        starting_line::IntOrNothing # block starting line index
        ending_line::IntOrNothing # block ending line index
        xpoints_per_line::IntOrNothing # number of X points in a single line
        ypoints_per_line::IntOrNothing# number of Y points in a single line
        decoding::Type{T} where T<:Decoding# type of single string decoding
        """
        JDXblock()
JDXreader obj constructor, creates empty object with no data
    """
        JDXblock()=begin
            new("",
                OrderedDict{String,Union{String,Float64}}(),
                Vector{Float64}(),
                Vector{Float64}(),
                Vector{Bool}(),
                nothing , nothing, nothing, nothing,NoDecoding)
        end
    end
    """
    write_jdx_file(file_name,jdx::JDXblock; kwargs...)

Writes [`JDXblock`](@ref) object to file
"""
write_jdx_file(file_name,jdx::JDXblock; kwargs...) = write_jdx_file(file_name,jdx.x_data,
                                                                        jdx.y_data,
                                                                        get(jdx.data_headers,"XUNITS","1/CM"),
                                                                        get(jdx.data_headers,"YUNITS","A.U.");
                                                                        kwargs...)
   """
    write_jdx_file(jdx::JDXblock; kwargs...)

Writes [`JDXblock`](@ref) object to file
"""
write_jdx_file(jdx::JDXblock; kwargs...) = write_jdx_file(jdx.file_name,jdx; kwargs...)                                                                  
    """
    read_jdx_file(file_name::String;delimiter=" ")

    Reads JCAMP-DX file

Input arguments: 

file_name - full file name
(optional) - delimiter 

returns namedtuple (or a vector of namedtuples in the case of multiple blocks) 

with fields:

x - coordinate (wavelength, wavenumber or other)

y - data

headers - dictionary in "String => value" format with headers values 

is_violated - is the Vector of Bool, if `is_violated[ind]` is `true` than control x-value
which stay at the `ind`'th  `line of data` in the file differs from the corresponding 
value of generated x-coordinate array by more than $(X_VIOLATION_CRITERIA). The `line of data`
index are the indces from 1 to `data_lines_number` where 
`data_lines_number=total_lines_number - (header_lines_number +1)`
`header_lines_number` is the number of headers keys, `total_lines_number` - is the total number of
lines in file (+1 stays to count the last line of file, which contains `#end=`  ) 
lines in file , see also [`read!`](@ref) and [`JDXblock`](@ref)

"""
    function read_jdx_file(file_name::String;delimiter = nothing)
        @assert isfile(file_name) "Must be a file"
        jdx_blocks = count_blocks(file_name)
        if length(jdx_blocks) ==1
            return read!(jdx_blocks[1],delimiter=delimiter)
        else
            return read!.(jdx_blocks,delimiter=delimiter)
        end
    end
    """
    JDXblock(file_name::String)

Creates JDXreader object from full file name 
"""
    function JDXblock(file_name::String) # external constructor
        if !isfile(file_name)
            error("Input filename must be the fullfile name to the existing file")
        end
        jdx = JDXblock()
        jdx.file_name = file_name
        return jdx
    end
    function parseJDXheaders(jdx::JDXblock,headers::Vector{String})
        for head in headers
            name = strip(head,'#')
            name = strip(name,'$')
            splitted_string = split(name,"=")
            splits_number = length(splitted_string)
            if splits_number==1
                jdx.data_headers[string(name)] = string(name)
                continue
            elseif splits_number==2
                jdx.data_headers[String(splitted_string[1])] = 
                isnothing(tryparse(Float64,splitted_string[2])) ? string(splitted_string[2]) : parse(Float64,splitted_string[2])  
            else
                jdx.data_headers[String(splitted_string[1])] = join(splitted_string[2:end],"=")
            end 
        end
    end

DEFAULT_DELIMITER(::Type{XYYline}) = isspace
DEFAULT_DELIMITER(::Type{XYXYline}) = r"[,;]" 

addline!(::Type{T},jdx::JDXblock, current_line;delimiter=DEFAULT_DELIMITER(T))  where T<:DATAline = addline!(T,
                                        jdx, 
                                        Float64[],
                                        current_line,
                                        1,delimiter=delimiter)
    """
    addline!(::Type{XYYline},jdx::JDXblock, data_chunk::AbstractVector,
                                current_line::String,
                                line_index::Int; # index of current data chunk
                                delimiter=isspace)

This function parses `current_line` string of file and fills the parsed data to the y-vector of `jdx`  object
`number_of_y_point_per_chunk` - the number of data point (excluding the x-coordinate) in the line 
`chunk_index`  - the index of chunk 
`delimiter`   - data points delimiter used in `split` function

"""
function addline!(::Type{XYYline},jdx::JDXblock, data_chunk::AbstractVector,
                                current_line::String,
                                line_index::Int; # index of current data chunk
                                delimiter=isspace) 
        
        if isempty(current_line) || current_line[1]=='#'  
             return NaN
        end 
        is_first_line = isnothing(jdx.ypoints_per_line)
        cur_points_number =  fill_data_chunk!(data_chunk,current_line,delimiter)
        if is_first_line
            number_of_y_point_per_chunk = cur_points_number - 1 
            jdx.ypoints_per_line = number_of_y_point_per_chunk
            jdx.xpoints_per_line = 1
            resize!(jdx.y_data,number_of_y_point_per_chunk)
        else
            number_of_y_point_per_chunk = jdx.ypoints_per_line    
        end 
        starting_index = 1 + (line_index-1) * number_of_y_point_per_chunk
        ending_index =  starting_index + cur_points_number - 2
        v = @view jdx.y_data[ starting_index : ending_index ]
        y_data_chunk =@view data_chunk[2:cur_points_number]
        copyto!(v,y_data_chunk)
        return data_chunk[1] # returns x-value for checks
    end
    function addline!(::Type{XYXYline},jdx::JDXblock, data_chunk::AbstractVector,
        current_line::String,
        line_index::Int; # index of current data chunk
        delimiter=r"[,;]") 

        if isempty(current_line) || current_line[1]=='#'  
                    return NaN
        end 
        is_first_line = isnothing(jdx.ypoints_per_line)
        cur_points_number =  fill_data_chunk!(data_chunk,current_line,delimiter)
        if is_first_line
            number_of_y_point_per_chunk = div(cur_points_number,2) 
            jdx.ypoints_per_line = number_of_y_point_per_chunk
            jdx.xpoints_per_line = cur_points_number - number_of_y_point_per_chunk
            resize!(jdx.y_data,number_of_y_point_per_chunk)
            resize!(jdx.x_data,jdx.xpoints_per_line)
        else
            number_of_y_point_per_chunk = jdx.ypoints_per_line    
        end 
        starting_index = 1 + (line_index-1) * number_of_y_point_per_chunk
        ending_index =  starting_index + div(cur_points_number,2) - 1
        vY = @view jdx.y_data[ starting_index : ending_index ]
        vX = @view jdx.x_data[ starting_index : ending_index ]
        y_data_chunk = @view data_chunk[2:2:cur_points_number]
        x_data_chunk = @view data_chunk[1:2:cur_points_number]
        copyto!(vY,y_data_chunk)
        copyto!(vX,x_data_chunk)
        return data_chunk[1] # returns x-value for checks
    end   
    """
    fill_data_chunk!(data_chunk::MVector,current_line,delimiter,chunk_counter::Int=1)

Fills vector from string line splitting it by the delimiter and parsing each value to Float64
Returns the number for parsed numbers
"""
function fill_data_chunk!(data_chunk::MVector,current_line,delimiter,chunk_counter::Int=1)
        #@show jdx.y_data
        for s in eachsplit(current_line,delimiter)
            s = strip(s)
            !isempty(s) || continue # check for empty string
            if !is_PAC_string(s) 
                data_chunk[chunk_counter] = Base.parse(Float64,s)
                chunk_counter += 1
            else
                out = split_PAC_string!(data_chunk,chunk_counter,s)
                chunk_counter += out
            end
        end
        return chunk_counter-1 # total number of parsed numerics
    end
    """
    fill_data_chunk!(data_chunk::Vector{Float64},current_line,delimiter,chunk_counter::Int=1)

Appends data to vector, initial chunk can be of zero size
"""
function fill_data_chunk!(data_chunk::Vector{Float64},current_line,delimiter,chunk_counter::Int=1)
        #@show jdx.y_data
        for s in eachsplit(current_line,delimiter)
            s = strip(s)
            !isempty(s) || continue # check for empty string
            if !is_PAC_string(s) 
                push!(data_chunk, Base.parse(Float64,s))
                chunk_counter += 1
            else
                parsed_numbers_vector = split_PAC_string(s)
                out = length(parsed_numbers_vector)
                append!(data_chunk,parsed_numbers_vector)
                chunk_counter += out
            end
        end
        return chunk_counter-1
    end
    """
    split_PAC_string(s::String, pattern::Regex=r"[+-]")

Function splits string with digits separated by multiple patterns and fills array `a`
starting from `starting_index` in-place, returns the quantity of numbers parsed from the string
"""
function split_PAC_string!(a::AbstractArray,starting_index::Int,s::AbstractString,
            pattern::Regex=r"[+-]")
        s = strip(s)
        offsets = [x.offset for x in eachmatch(pattern, s)] # generating array of patterns match
        N = length(offsets)
        if N==0 # there are no patterns
            val =  Base.tryparse(Float64,s)
            !isnothing(val) ? a[1]=val : return starting_index+1
            return starting_index+1
        end
        if offsets[1]==1
            starting_counter  = 2
        else
            starting_counter = 1
        end
        counter=starting_index
        start_index = 1
        for ii in starting_counter:N
            stop_index = offsets[ii]-1
            a[counter] = Base.parse(Float64,s[start_index:stop_index]) 
            start_index = 1 + stop_index 
            counter +=1
        end
        a[counter] = Base.parse(Float64,s[start_index:end])
        return counter-starting_index+1
    end
    """
    split_PAC_string(s::String, pattern::Regex=r"[+-]")

Function splits string with digits separated by multiple patterns
"""
function split_PAC_string(s::AbstractString, pattern::Regex=r"[+-]")
        s = strip(s)
        offsets = [x.offset for x in eachmatch(pattern, s)]
        N = length(offsets)
        parts = Float64[]
        if N==0 
            val =  Base.tryparse(Float64,s)
            !isnothing(val) ? push!(parts, Base.parse(Float64,s)) : return parts
            return parts
        end
        if offsets[1]==1
            starting_counter  = 2
        else
            starting_counter = 1
        end
        start_index = 1
        for ii in starting_counter:N
            stop_index = offsets[ii]-1
            push!(parts, Base.parse(Float64,s[start_index:stop_index])) 
            start_index = 1 + stop_index 
        end
        push!(parts, Base.parse(Float64,s[start_index:end])) 
        return parts
    end

    """
    generateXvector!(jdx::JDXblock)

Generates equally spaced x-vector 
"""
function generateXvector!(jdx::JDXblock)
            point_number = round(Int64,jdx.data_headers["NPOINTS"])
            starting_X = haskey(jdx.data_headers,"FIRSTX") ? jdx.data_headers["FIRSTX"] : 0.0
            if haskey(jdx.data_headers,"XFACTOR")
                step_value =  ((jdx.data_headers["LASTX"] -  starting_X)/jdx.data_headers["XFACTOR"])/(point_number-1)
            else
                step_value =  (jdx.data_headers["LASTX"] -  starting_X)/(point_number-1)
            end
            if haskey(jdx.data_headers,"DELTAX") # check for the delta x values
                delta_step = abs(1 - jdx.data_headers["DELTAX"]/step_value)
                delta_step < 1e-3 ? nothing : @warn "DELTAX value for TITLE=$(jdx.data_headers["TITLE"]) is violated (DELTAX-ACTUAL_STEP)/ACTUAL_STEP = $(delta_step) should be less than $(DELTAX_VIOLATION_CRITERIUM)"
            end
            resize!(jdx.x_data,point_number)
            map!(i->starting_X + i*step_value,jdx.x_data,0:point_number-1)
    end

    """
    parse_headers(file::String)

Parses headers from JCAMP-DX file
"""
parse_headers(file::String) = file |> JDXblock |>  parse_headers!

    function parse_headers!(jdx::JDXblock)
        read!(jdx; only_headers=true)
        return jdx.data_headers
    end

    function is_supported_jdx_format(file_name::AbstractString)
        headers = parse_headers(file_name)
        return haskey(headers,"JCAMP-DX") && isa(headers["JCAMP-DX"],Number) && headers["JCAMP-DX"]<=SUPPORTED_JCAMPDX_VERSION
    end

    """
    read!(jdx::JDXblock)

fills precreated JDXblock object see [`JDXblock`](@ref)
"""
function read!(jdx::JDXblock; delimiter=nothing,
                              only_headers::Bool=false)
        if !isfile(jdx.file_name)
            return nothing
        end
        header_lines = Vector{String}()
        isnothing(jdx.starting_line) || jdx.starting_line<0 ? jdx.starting_line = 1 : nothing
        isnothing(jdx.ending_line) || jdx.ending_line < 0 ? jdx.ending_line = countlines(jdx.file_name) : nothing 
        total_number_of_lines = jdx.ending_line - jdx.starting_line + 1
        x_point =0.0;
        line_type =  nothing
        #jdx.y_data = Vector{Float64}()
        open(jdx.file_name) do io_file
            x_point = 0.0
            is_data_started = false # tru if data block is already started
            header_lines_counter = 0
            
            if !isnothing(jdx.starting_line) && jdx.starting_line>1
                skip_lines = jdx.starting_line - 1
            else
                skip_lines = 0
            end
            for ln in Iterators.drop(eachline(io_file), skip_lines)#eachline(io_file) #scans only headers and the first line of data
                if ln[1]!='#' 
                    if is_data_started 
                        #  reading the first line to get the number of data points per line et.c.
                        if !only_headers
                            delimiter =  isnothing(delimiter) ? DEFAULT_DELIMITER(line_type) : delimiter
                            jdx.decoding = get_decoding_type(ln) # filling decoding type
                            x_point = addline!(line_type, jdx, ln, delimiter = delimiter) # this function modifies
                            # jdx file block by setting values to the number of points-per-line properties 
                        end
                        break 
                        #  parsing the first line
                    else #  it is possible that the header is multilined, in this case the data is not started, but line has no leading "#" 
                        # in this case, multilined header is parsed to a single line 
                        ln_total = header_lines[end] # taking previous vector element as a starting content of this line
                        ln_cur=""
                        for ln_cur in eachline(io_file)
                            ln_cur[1]!='#' ?  ln_total *=" "*ln_cur : break
                            header_lines_counter+=1
                        end
                        push!(header_lines,ln_total)
                        header_lines_counter+=1
                        push!(header_lines,ln_cur)
                    end
                else
                    if !is_data_started 
                        l_upped = uppercase(ln)
                        is_data_started = occursin("##XYDATA",l_upped) || occursin("##XYPOINTS",l_upped) || occursin("##PEAK TABLE",l_upped)
                        if is_data_started
                            ln = l_upped
                            if occursin("(X++(Y..Y))",ln)
                                line_type =  XYYline
                            elseif  occursin("(XY..XY)",ln)
                                line_type =  XYXYline
                            end
                        end
                    end 
                    # check if data is already started
                    header_lines_counter+=1
                    push!(header_lines,ln)
                end
            end # endof header lines parsing 

            parseJDXheaders(jdx,header_lines) #perses headers to dict from a vector of strings
            if only_headers # if we need only block headers without data 
                return (x=Vector{Float64}([]),
                        y=Vector{Float64}([]), 
                        headers=jdx.data_headers,
                        is_violated = false)
            end
            number_of_y_point_per_chunk = length(jdx.y_data) # number of numeric points per line
            # by default we assume that the last line is ##END= thus it is ignored
            data_lines_number = total_number_of_lines - header_lines_counter - 1 # last line is ignored
            
            if haskey(jdx.data_headers,"NPOINTS") # correct JDX file
                total_point_number = round(Int64,jdx.data_headers["NPOINTS"]) # total number of points is known
                jdx.is_violated_flag = Vector{Bool}(undef,data_lines_number)
                is_XYYline = line_type <: XYYline
                resize!(jdx.y_data,total_point_number)
                is_XYYline ? generateXvector!(jdx) : resize!(jdx.x_data,total_point_number)
                points_number_per_chunk = jdx.ypoints_per_line + jdx.xpoints_per_line
                data_chunk_container = MVector{points_number_per_chunk,Float64}(undef)    
                if is_XYYline    
                    x_gen = jdx.x_data[1] # generated x to compare
                    jdx.is_violated_flag[1] = !isapprox(x_gen, x_point,rtol=X_VIOLATION_CRITERIA)    
                    for i in 2:data_lines_number
                        ln = jdx.decoding(readline(io_file))
                        x_point = addline!(line_type,jdx,data_chunk_container,ln,i,delimiter=delimiter)
                        if !isnan(x_point)
                            x_gen = jdx.x_data[(i-1)*number_of_y_point_per_chunk + 1] # generated x
                            jdx.is_violated_flag[i] =  !isapprox(x_gen, x_point,rtol=X_VIOLATION_CRITERIA)
                        else
                            l = length(jdx.is_violated_flag)
                            deleteat!(jdx.is_violated_flag,i:l)
                            break
                        end#endif       
                    end#endfor
                else #XYXYline
                   #@show jdx
                    for i in 2:data_lines_number
                        ln = jdx.decoding(readline(io_file))
                        addline!(line_type,jdx,data_chunk_container,ln,i,delimiter=delimiter)     
                    end
                end
            else number_of_y_point_per_chunk==1 # file with two columns like CSV
                total_point_number = data_lines_number
                jdx.is_violated_flag = Vector{Bool}(undef,1)
                resize!(jdx.y_data,total_point_number)
                jdx.x_data = similar(jdx.y_data)
                jdx.x_data[1] = x_point
                data_chunk = MVector{2,Float64}(undef)
                for i in 2:data_lines_number
                    ln = readline(io_file)
                     x_out = addline!(XYXYline,jdx,data_chunk, ln,i,delimiter=delimiter)
                     !isnan(x_out) ? jdx.x_data[i] = x_out : nothing
                end
                jdx.is_violated_flag[1]=false
            end

        end # close file
        if line_type<:XYYline   
            is_violated = any(jdx.is_violated_flag)
            !is_violated ? nothing : @warn "X values check for TITLE=$(jdx.data_headers["TITLE"]) is violated at $(sum(jdx.is_violated_flag)) points"
        end
        if haskey(jdx.data_headers,"YFACTOR") 
            y_factor::Float64 = jdx.data_headers["YFACTOR"] 
            jdx.y_data .*= y_factor
        end        
        return (x=jdx.x_data,
                y=jdx.y_data, 
                headers=jdx.data_headers,
                is_violated = jdx.is_violated_flag)
    end#enof read!

    function count_blocks(file_name)
        @assert isfile(file_name) "Not a file"
        jdx_blocks = Vector{JDXblock}()
        open(file_name) do io
            block_started = false
            for (line_counter,ln) in enumerate(eachline(io))
                #@show line_counter,ln
                if ln[1] =='#'
                    if (!block_started) && occursin("##TITLE",uppercase(ln))
                        jdx_block = JDXblock(file_name)
                        jdx_block.starting_line = line_counter
                        push!(jdx_blocks,jdx_block)
                        block_started = true    
                    end
                    if block_started && occursin("##END",uppercase(ln))
                        jdx_blocks[end].ending_line = line_counter
                        block_started = false    
                    end
                end        
            end#forloop

        end # fileopen
        return jdx_blocks
    end#funend

end
