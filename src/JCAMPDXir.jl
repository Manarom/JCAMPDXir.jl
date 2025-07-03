
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
    const IntOrNothing = Union{Int,Nothing}
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
    include("DataViolation.jl")
    include("WriteJCAMP.jl")
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
By default there is no decoding `No_Decoding`, but if the first line of data 
in block contains `SQZ_digits` than it marks entire block as `SQZ_Decoding`
by setting the decoding field of `JDXblock` object. String decoding is performed 
right after the reading the line of data from file.

    """
    abstract type Decoding end
    struct SQZ_Decoding<:Decoding end
    struct No_Decoding<:Decoding end
    struct DIF_Decoding<:Decoding end
    struct DUP_Decoding<:Decoding end
    struct DIF_DUP_Decoding<:Decoding end
    struct SQZ_DIF_Decoding<:Decoding end
    struct SQZ_DUP_Decoding<:Decoding end
    struct SQZ_DIF_DUP_Decoding<:Decoding end
    const ASDF = Union{SQZ_DUP_Decoding,SQZ_DIF_DUP_Decoding,DUP_Decoding}
# this functions should return the string with all values replaced by the actual numers
    decode(::Type{No_Decoding},s::AbstractString) = s
    decode(::Type{SQZ_Decoding},s) = replace(s,SQZ_digits...)
    decode(::Type{DIF_Decoding},s) = begin 
        s = replace(s,DIF_digits...)
    end    
    #const ASDF = Union{}

    decode(::Type{T},s) where T<:ASDF= begin
        s = decode(SQZ_Decoding,s) # decoding 
        return decode(DIF_Decoding,s)
    end
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

"""
    get_decoding_type(s::AbstractString)

Returns decoding type by parsing the input string
"""
function get_decoding_type(s::AbstractString)
    is_SQZ = is_SQZ_string(s)
    is_DIF = is_DIF_string(s)
    is_DUP = is_DUP_string(s)
    is_SQZ && is_DIF && is_DUP && return SQZ_DIF_DUP_Decoding
    is_SQZ && is_DIF && return SQZ_DIF_Decoding
    is_SQZ && is_DUP && return SQZ_DUP_Decoding
    is_DIF && is_DUP && return DIF_DUP_Decoding
    is_SQZ && return SQZ_Decoding
    is_DIF && return DIF_Decoding
    is_DUP && return DUP_Decoding

    return No_Decoding
end

    function Base.parse(::Type{T},s::AbstractString) where T<:Union{SQZ_DIF_DUP_Decoding,DUP_Decoding,DIF_DUP_Decoding}
        is_DUP_string(s) || return [Base.parse(Float64,s)]
        m = match(DUP_regexp,s)
        !isnothing(m) || return Float64[]
        n = DUP_digits[m.match[1]] # multiplyer numeric
        return  Base.parse(Float64,s[1:m.offset-1])
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
    Stored parsed data
Must be filled using [`read!`](@ref) function
"""
    mutable struct JDXblock
        # Main struct, prepares  file name, parses data and data headers
        file_name::String
        data_headers::Dict{String,Union{String,Float64}}
        x_data::Vector{Float64}
        y_data::Vector{Float64}
        violation::DataViolation# 
        starting_line::IntOrNothing # block starting line index
        ending_line::IntOrNothing # block ending line index
        #xpoints_per_line::IntOrNothing # number of X points in a single line
        #ypoints_per_line::IntOrNothing# number of Y points in a single line
        filled_points_counter::Int # stores the number of filled points 
        #decoding# decoding function 
        """
        JDXblock()
JDXreader obj constructor, creates empty object with no data
    """
        JDXblock()=begin
            new("",#filename
                OrderedDict{String,Union{String,Float64}}(),#headers
                Vector{Float64}(),# x
                Vector{Float64}(),# y
                DataViolation(),#is_violated
                nothing ,# starting_line
                nothing, # ending_line
                0)#,
                #No_Decoding)
        end
    end
    # 
    mutable struct DataBuffer{  DataLineType<:DATAline, # line format XYYline or XYXYline
                                BufferType<:AbstractVector,# data buffer type may be static vector or dynamic vector
                                LineDecodingType<:Decoding, # this type performs some work on line before parsing 
                                ParseToType<:Union{Float64,ASDF} # this is a type of single chunk to be prased to using Base.parse
                                }
        buffer::BufferType
        xpoints_per_line::IntOrNothing # number of X points in a single line
        ypoints_per_line::IntOrNothing# number of Y points in a single line
        DataBuffer() = new{XYYline,Vector{Float64},No_Decoding,Float64}(Float64[],nothing,nothing)
        function DataBuffer(buffer::BufferType;
                            DataLineType = XYYline,
                            LineDecodingType = No_Decoding,
                            ParseToType = Float64,
                            xpoints_per_line::IntOrNothing=nothing,
                            ypoints_per_line::IntOrNothing=nothing) where BufferType<:AbstractVector 

            return new{DataLineType,BufferType,LineDecodingType,ParseToType}(buffer,xpoints_per_line,ypoints_per_line)
        end
    end
    is_resizable(db::DataBuffer) = buffertype(db) <:Vector
    datalinetype(::DataBuffer{DataLineType}) where DataLineType = DataLineType
    buffertype(::DataBuffer{DL,BufferType}) where {DL,BufferType} = BufferType
    linedecodingtype(::DataBuffer{DL,BT,LineDecodingType}) where {DL,BT,LineDecodingType} = LineDecodingType
    parsetotype(::DataBuffer{DL,BT,LD,ParseToType}) where {DL,BT,LD,ParseToType} = ParseToType
    Base.length(db::DataBuffer) = length(db.buffer)
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
value of generated x-coordinate array by more than $(X_VIOLATION_CRITERIUM). The `line of data`
index are the indces from 1 to `data_lines_number` where 
`data_lines_number=total_lines_number - (header_lines_number +1)`
`header_lines_number` is the number of headers keys, `total_lines_number` - is the total number of
lines in file (+1 stays to count the last line of file, which contains `#end=`  ) 
lines in file , see also [`read!`](@ref) and [`JDXblock`](@ref)

"""
    function read_jdx_file(file_name::String;fixed_columns_number::Bool=false,delimiter = nothing)
        @assert isfile(file_name) "Must be a file"
        jdx_blocks = count_blocks(file_name)
        if length(jdx_blocks) ==1
            return read!(jdx_blocks[1],delimiter=delimiter,fixed_columns_number=fixed_columns_number)
        else
            return read!.(jdx_blocks,delimiter=delimiter,fixed_columns_number=fixed_columns_number)
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
        multilined_header_started = false
        last_key = ""
        for head in headers
            !isempty(head) || continue
            if head[1]=='#'
                name = strip(head,'#')
                name = strip(name,'$')
                splitted_string = split(name,"=")
                splits_number = length(splitted_string)
                if splits_number==1
                    last_key = string(name)
                    last_val = string(name)
                elseif splits_number==2
                    last_key = String(splitted_string[1])
                    last_val = isnothing(tryparse(Float64,splitted_string[2])) ? string(splitted_string[2]) : parse(Float64,splitted_string[2]) 
                else
                    last_key = String(splitted_string[1])    
                    last_val = join(splitted_string[2:end],"=")
                end 
                jdx.data_headers[last_key] = last_val
            else
                val = String(jdx.data_headers[last_key])
                jdx.data_headers[last_key] = val *" \n "*head
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
    addline!(::Type{XYYline},jdx::JDXblock, data_buffer::AbstractVector,
                                current_line::String,
                                line_index::Int; # index of current data chunk
                                delimiter=isspace)

This function parses `current_line` string of file and fills the parsed data to the y-vector of `jdx`  object
`number_of_y_point_per_chunk` - the number of data point (excluding the x-coordinate) in the line 
`chunk_index`  - the index of chunk 
`delimiter`   - data points delimiter used in `split` function

"""
    function addline!( jdx::JDXblock, 
                        data_buffer::DataBuffer{DataLineType},
                        current_line::String,
                        line_index::Int; # index of current data chunk
                        delimiter=isspace) where DataLineType<:XYYline
        
        if isempty(current_line) || current_line[1]=='#'  
             return NaN
        end 
        is_first_line = jdx.filled_points_counter == 0
        cur_points_number =  fill_data_buffer!(data_buffer,current_line,delimiter)
        if is_first_line
            number_of_y_point_per_chunk = cur_points_number - 1 
            data_buffer.ypoints_per_line = number_of_y_point_per_chunk
            data_buffer.xpoints_per_line = 1
            #jdx.ypoints_per_line = number_of_y_point_per_chunk
            #jdx.xpoints_per_line = 1
            resize!(jdx.y_data,number_of_y_point_per_chunk)
        end 
        starting_index = 1 + jdx.filled_points_counter
        ending_index =  starting_index + cur_points_number - 2
        v = @view jdx.y_data[ starting_index : ending_index ]
        y_data_buffer =@view data_buffer.buffer[2:cur_points_number]
        copyto!(v,y_data_buffer)
        jdx.filled_points_counter += cur_points_number - 1 
        return data_buffer.buffer[1] # returns x-value for checks
    end
    function addline!(::Type{XYXYline},
        jdx::JDXblock, 
        data_buffer::AbstractVector,
        current_line::String,
        line_index::Int; # index of current data chunk
        delimiter=r"[,;]") 

        if isempty(current_line) || current_line[1]=='#'  
                    return NaN
        end 
        is_first_line = jdx.filled_points_counter == 0
        cur_points_number =  fill_data_buffer!(data_buffer,current_line,delimiter)
        if is_first_line
            number_of_y_point_per_chunk = div(cur_points_number,2) 
            jdx.ypoints_per_line = number_of_y_point_per_chunk
            jdx.xpoints_per_line = cur_points_number - number_of_y_point_per_chunk
            resize!(jdx.y_data,number_of_y_point_per_chunk)
            resize!(jdx.x_data,jdx.xpoints_per_line)   
        end 
        starting_index = 1 + jdx.filled_points_counter
        ending_index =  starting_index + div(cur_points_number,2) - 1
        vY = @view jdx.y_data[ starting_index : ending_index ]
        vX = @view jdx.x_data[ starting_index : ending_index ]
        y_data_buffer = @view data_buffer[2:2:cur_points_number]
        x_data_buffer = @view data_buffer[1:2:cur_points_number]
        copyto!(vY,y_data_buffer)
        copyto!(vX,x_data_buffer)
        jdx.filled_points_counter += div(cur_points_number,2) 
        return data_buffer[1] # returns x-value for checks
    end   
    """
    fill_data_buffer!(data_buffer::MVector,current_line,delimiter,chunk_counter::Int=1)

Fills vector from string line splitting it by the delimiter and parsing each value to Float64
Returns the number for parsed numbers
"""
function fill_data_buffer!(data_buffer::DataBuffer{DataLineType,BufferType,LineDecodingType,ParseToType},
                                                current_line,
                                                delimiter,
                                                chunk_counter::Int=1) where {DataLineType, 
                                                                BufferType<:MVector,
                                                                LineDecodingType,
                                                                ParseToType}
                                                # DataBuffer{BufferType,ParseToType,DataLineType,LineDecodingType}
        current_line = LineDecodingType(current_line)
        for s in eachsplit(current_line,delimiter)
            s = strip(s)
            !isempty(s) || continue # check for empty string
            if chunk_counter>length(data_buffer)
                break
            end
            if !is_PAC_string(s) && ParseToType<:Float64
                data_buffer.buffer[chunk_counter] = Base.parse(ParseToType,s)
                chunk_counter += 1
            else
                out = split_PAC_string!(data_buffer,chunk_counter,s)
                chunk_counter += out
            end
        end
        return chunk_counter-1 # total number of parsed numerics
    end
    """
    fill_data_buffer!(data_buffer::Vector{Float64},current_line,delimiter,chunk_counter::Int=1)

Appends data to vector, initial chunk can be of zero size

"""
function fill_data_buffer!(data_buffer::DataBuffer{DataLineType,BufferType,LineDecodingType,ParseToType},
                                            current_line,
                                            delimiter,
                                            chunk_counter::Int=1) where {DataLineType, 
                                            BufferType<:Vector,
                                            LineDecodingType,
                                            ParseToType}
        #@show jdx.y_data
        current_line = LineDecodingType(current_line)
        for s in eachsplit(current_line,delimiter)
            s = strip(s)
            !isempty(s) || continue # check for empty string
            if !is_PAC_string(s) 
                if ParseToType<:Float64
                    set_or_push!(data_buffer.buffer,chunk_counter,Base.parse(ParseToType,s))
                    chunk_counter += 1
                else
                    error("Under construction")
                end
            else
                out = split_PAC_string!(data_buffer,chunk_counter,s)
                chunk_counter += out
            end
        end
        points_in_chunk = chunk_counter-1
        length(data_buffer)!=points_in_chunk &&  resize!(data_buffer.buffer,points_in_chunk)
        return points_in_chunk
    end

    function set_or_push!(v::Vector{T}, i::Int, val::T) where T<:Number
        if 1 <= i <= length(v)
            v[i] = val
        elseif i == length(v) + 1
            push!(v, val)
        else i > length(v) + 1
            resize!(v, i)  # resize to accommodate the index
            v[i] = val
        end
    end
    function set_or_push!(dest::Vector{T},i::Int,val::AbstractVector{T}) where T<:Number
        end_index = i + length(val)-1
        for ii in end_index:-1:i
            v_ind =  ii -i + 1
            set_or_push!(dest,ii,val[v_ind])
        end
    end
"""
    split_PAC_string!(a::Vector,starting_index::Int,s::AbstractString,
        pattern::Regex=r"[+-]")

Version of splitting PCA string which fills resizable vector a
"""
function split_PAC_string!(a::Vector{T},starting_index::Int,s::AbstractString,
        pattern::Regex=r"[+-]") where T<:Number
    starting_index > length(a) && resize!(a,starting_index)
    
    s = strip(s)
    counter = starting_index
    stop_index = 0
    stop_index = 0
    start_index = 1
    for (ii,m) in enumerate(eachmatch(pattern,s))
        offset = getfield(m,:offset)
        ii==1 && offset == 1 && continue
        stop_index = offset-1
        set_or_push!(a,counter, Base.parse(T,s[start_index:stop_index])) 
        start_index = 1 + stop_index 
        counter +=1
    end   
    set_or_push!(a,counter, Base.parse(T,s[start_index:end]))
    return counter-starting_index+1
end
    """
    split_PAC_string(s::String, pattern::Regex=r"[+-]")

Function splits string with digits separated by multiple patterns and fills array `a`
starting from `starting_index` in-place, returns the quantity of numbers parsed from the string
"""
function split_PAC_string!(a::MVector{N,T},starting_index::Int,s::AbstractString,
            pattern::Regex=r"[+-]") where {N,T}
        s = strip(s)
        counter = starting_index
        stop_index = 0
        stop_index = 0
        start_index = 1
        for (ii,m) in enumerate(eachmatch(pattern,s))
            offset = getfield(m,:offset)
            ii==1 && offset == 1 && continue
            stop_index = offset-1
            a[counter] = Base.parse(T,s[start_index:stop_index])
            start_index = 1 + stop_index 
            counter +=1
        end   
        a[counter] =  Base.parse(T,s[start_index:end])
        return counter-starting_index+1
    end
    """
    generateXvector!(jdx::JDXblock)

Generates equally spaced x-vector 
"""
function generateVectors!(jdx::JDXblock,::Type{XYYline})
            point_number = round(Int64,jdx.data_headers["NPOINTS"])
            starting_X = haskey(jdx.data_headers,"FIRSTX") ? jdx.data_headers["FIRSTX"] : 0.0
            if haskey(jdx.data_headers,"XFACTOR")
                step_value =  ((jdx.data_headers["LASTX"] -  starting_X)/jdx.data_headers["XFACTOR"])/(point_number-1)
            else
                step_value =  (jdx.data_headers["LASTX"] -  starting_X)/(point_number-1)
            end
            if haskey(jdx.data_headers,"DELTAX") # check for the delta x values
                check_data_point!(jdx.violation,jdx.data_headers["DELTAX"],
                                    step_value,:deltax)
            end
            resize!(jdx.x_data,point_number)
            map!(i->starting_X + i * step_value, jdx.x_data,0 : point_number-1)
            resize!(jdx.y_data,point_number)
    end
function generateVectors!(jdx::JDXblock,::Type{XYXYline})
        point_number = round(Int64,jdx.data_headers["NPOINTS"])
        resize!(jdx.x_data,point_number)
        resize!(jdx.y_data,point_number)
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
    read!(jdx::JDXblock; delimiter=nothing,
                                    only_headers::Bool=false,
                                    fixed_numbers_line::Bool=true)

fills precreated JDXblock object see [`JDXblock`](@ref)
"""
function read!(jdx::JDXblock; delimiter=nothing,
                               only_headers::Bool=false,
                               fixed_columns_number::Bool=true,
                               fixed_data_decoding::Bool = true)
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
            if !isnothing(jdx.starting_line) && jdx.starting_line>1
                skip_lines = jdx.starting_line - 1
            else
                skip_lines = 0
            end
            is_data_started = false # tru if data block is already started
            header_lines_counter = 0 # conter of lines per header
            first_line_buffer = DataBuffer()
            for ln in Iterators.drop(eachline(io_file), skip_lines)#eachline(io_file) #scans only headers and the first line of data
                # !isempty(ln) || continue
                if is_data_started 
                        if !only_headers
                            delimiter =  isnothing(delimiter) ? DEFAULT_DELIMITER(line_type) : delimiter
                            #jdx.decoding =  # filling decoding type
                            # @show jdx.decoding 

                            first_line_buffer = DataBuffer(Vector{Float64}(undef,6),
                                                        DataLineType =line_type,
                                                        LineDecodingType = get_decoding_type(ln) )
                            x_point = addline!( jdx, first_line_buffer,ln,1, delimiter = delimiter) # this function modifies
                            # jdx file block by setting values to the number of points-per-line properties 
                        end
                        break 
                else
                    l_upped = uppercase(ln)
                    is_data_started = occursin("##XYDATA",l_upped) || occursin("##XYPOINTS",l_upped) || occursin("##PEAK TABLE",l_upped)
                    if is_data_started
                        if occursin("(X++(Y..Y))",l_upped)
                            line_type =  XYYline
                        elseif  occursin("(XY..XY)",l_upped)
                            line_type =  XYXYline
                        else
                            error("Unknown line type")
                        end
                    end
                    # check if data is already started
                    header_lines_counter+=1
                    push!(header_lines,l_upped)
                end
            end # endof header lines parsing 
           # #@show header_lines
            parseJDXheaders(jdx,header_lines) #parses headers to dict from a vector of strings
            if only_headers # if we need only block headers without data 
                return (x=Vector{Float64}([]),
                        y=Vector{Float64}([]), 
                        headers=jdx.data_headers,
                        is_violated = false)
            end
            number_of_y_point_per_chunk = length(jdx.y_data) # number of numeric points per line
            # by default we assume that the last line is ##END= thus it is ignored
            data_lines_number = total_number_of_lines - header_lines_counter - 1 # last line is ignored because there should be ##END
            # @show jdx.data_headers["NPOINTS"]
            # @show jdx.ypoints_per_line
            if haskey(jdx.data_headers,"NPOINTS") # correct JDX file
                total_point_number = round(Int64,jdx.data_headers["NPOINTS"]) # total number of points is known
                is_XYYline = line_type <: XYYline
                #resize!(jdx.y_data,total_point_number)
                generateVectors!(jdx,line_type)
                #is_XYYline ? generateXvector!(jdx) : resize!(jdx.x_data,total_point_number)
                points_number_per_chunk = length(first_line_buffer)
                #@show points_number_per_chunk
                if !fixed_columns_number
                    data_buffer = first_line_buffer
                else
                    data_buffer = DataBuffer(MVector{points_number_per_chunk,Float64}(undef),
                                             DataLineType=line_type,
                                             LineDecodingType = linedecodingtype(first_line_buffer),
                                             ParseToType = parsetotype(first_line_buffer),
                                             xpoints_per_line = first_line_buffer.xpoints_per_line,
                                             ypoints_per_line = first_line_buffer.ypoints_per_line)
                end    
                # @show data_buffer_container
                x_factor = get(jdx.data_headers,"XFACTOR",1.0)
                if is_XYYline    
                    x_gen = jdx.x_data[1] # generated x to compare
                    #jdx.is_violated_flag[1] = !isapprox(x_gen, x_factor*x_point,rtol=X_VIOLATION_CRITERIUM)    
                    for i in 2:data_lines_number
                        ln = readline(io_file)
                        # @show length(jdx.y_data)
                        x_point = addline!(jdx,data_buffer,ln,i,delimiter=delimiter)
                        x_point = x_point * x_factor
                        if !isnan(x_point) #&& fixed_columns_number
                            x_gen = jdx.x_data[(i-1)*number_of_y_point_per_chunk + 1] # generated x
                            #jdx.is_violated_flag[i] =  !isapprox(x_gen, x_point,rtol=X_VIOLATION_CRITERIUM)
                        else
                            l = length(jdx.is_violated_flag)
                            #deleteat!(jdx.is_violated_flag,i:l)
                            break
                        end#endif       
                    end#endfor
                else #XYXYline
                   # @show jdx
                    for i in 2:data_lines_number
                        ln = jdx.decoding(readline(io_file))
                        addline!(line_type,jdx,data_buffer_container,ln,i,delimiter=delimiter)     
                    end
                end
            else number_of_y_point_per_chunk==1 # file with two columns like CSV
                total_point_number = data_lines_number
                jdx.is_violated_flag = Vector{Bool}(undef,1)
                resize!(jdx.y_data,total_point_number)
                jdx.x_data = similar(jdx.y_data)
                jdx.x_data[1] = x_point
                data_buffer = MVector{2,Float64}(undef)
                for i in 2:data_lines_number
                    ln = readline(io_file)
                     x_out = addline!(XYXYline,jdx,data_buffer, ln,i,delimiter=delimiter)
                     !isnan(x_out) ? jdx.x_data[i] = x_out : nothing
                end
                jdx.is_violated_flag[1]=false
            end

        end # close file
        if line_type<:XYYline   
            #is_violated = any(jdx.is_violated_flag)
           #!is_violated ? nothing : @warn "X values check for TITLE=$(jdx.data_headers["TITLE"]) is violated at $(sum(jdx.is_violated_flag)) points"
        end
        if haskey(jdx.data_headers,"YFACTOR") 
            y_factor::Float64 = jdx.data_headers["YFACTOR"] 
            jdx.y_data .*= y_factor
        end        
        return (x=jdx.x_data,
                y=jdx.y_data, 
                headers=jdx.data_headers)
    end#enof read!

    function count_blocks(file_name)
        @assert isfile(file_name) "Not a file"
        jdx_blocks = Vector{JDXblock}()
        open(file_name) do io
            block_started = false
            for (line_counter,ln) in enumerate(eachline(io))
                # @show line_counter,ln
                if !isempty(ln) && ln[1] =='#'
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
