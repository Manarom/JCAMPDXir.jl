
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
    """
    Squezzed symbols to string digits converter

    """
    const SQZ_digits = Dict('@'=>" 0", 'A'=>" 1", 'B'=>" 2", 
                            'C'=>" 3", 'D'=>" 4",'E'=>" 5",
                            'F'=>" 6", 'G'=>" 7", 'H'=>" 8", 
                            'I'=>" 9",'a'=>" -1", 'b'=>" -2", 
                            'c'=>" -3", 'd'=>" -4", 'e'=>" -5", 
                            'f'=>" -6", 'g'=>" -7", 'h'=>" -8", 
                            'i'=>" -9")#','=>" " '+'=>"+",  '-'=>"-",  
                            #','=>" ") misleading with XY...XY delimiter and PAC form
    """
    Differential symbols to string digits converter

    """
	const DIF_digits = Dict('%'=>"0", 'J'=>"1",  'K'=>"2",  
                            'L'=>"3",  'M'=>"4",  'N'=>"5",  
                            'O'=>"6",  'P'=>"7",  'Q'=>"8",
                            'R'=>"9",  'j'=>"-1", 'k'=>"-2", 
                            'l'=>"-3", 'm'=>"-4", 'n'=>"-5",
                            'o'=>"-6", 'p'=>"-7", 'q'=>"-8",
                            'r'=>"-9")

    const DIF_shift = Dict(k=>" "*string(k) for k in keys(DIF_digits))
    """
    Duplication symbols to string digits

    """
    const DUP_digits = Dict('S'=>1, 'T'=>2, 'U'=>3, 'V'=>4, 
                            'W'=>5, 'X'=>6, 'Y'=>7, 'Z'=>8,
                            's'=>9)
    const DUP_shift = Dict(k=>string(k)*" " for k in keys(DUP_digits))

    const SQZ_regexp = Regex("["*join(keys(SQZ_digits))*"]" )
    const DIF_regexp = Regex("["*join(keys(DIF_digits))*"]" )
    const DUP_regexp = Regex("["*join(keys(DUP_digits))*"]" )

    """
    Types for decoding of both lines and chunks, line and chunk typisation is used 
    to dispath during parsing
"""
    abstract type Decoding end
    """
    (::Type{T})(s::AbstractString) where T<:Decoding

By default, all chunks and line decoding do nothing to the string line 

"""
    (::Type{T})(s::AbstractString) where T<:Decoding= s   
    """
    Types for decoding lines and shunk string
"""
    abstract type LineDecoding<:Decoding end
    """
    When there is no line decoding all 
"""
    struct No_Line_Decoding<:LineDecoding end
    struct SQZ<:LineDecoding end
    struct PAC<:LineDecoding end
    struct SQZ_PAC<:LineDecoding end
    """
    When using Unspecified_Line type line type parser looks for the actual type for each line

    """
    struct Unspecified_Line<:LineDecoding end
# string type checkers
    is_PAC_string(s::AbstractString) =occursin('-',s)||occursin('+',s)
    is_SQZ_string(s::AbstractString) = occursin(SQZ_regexp,s)
    is_SQZ_PAC_string(s::AbstractString) = is_PAC_string(s) && is_SQZ_string(s)
    """
    get_line_decoding(s::AbstractString)

Returns line type by searching for specific symbols
"""
function get_line_decoding(s::AbstractString)
        is_PAC = is_PAC_string(s)
        is_SQZ = is_SQZ_string(s)
        is_PAC && is_SQZ && return SQZ_PAC
        is_PAC && return PAC
        is_SQZ && return SQZ
        return No_Line_Decoding
    end

    """
    (::Type{SQZ})(s::AbstractString)

All types are callable, when calling on a string 
"""
(::Type{SQZ})(s::AbstractString) = replace(s,SQZ_digits...)
    (::Type{PAC})(s::AbstractString) = replace(s,"+"=>" +","-"=>" -")
    (::Type{SQZ_PAC})(s::AbstractString) = s |> PAC |> SQZ
    (::Type{Unspecified_Line})(s::AbstractString) = get_line_decoding(s)(s)

# chunk decoding types
    abstract type ChunkDecoding <: Decoding end
    struct No_Chunk_Decoding<: ChunkDecoding end
    struct DIF <: ChunkDecoding end
    struct DUP <: ChunkDecoding end
    struct DIF_DUP <: ChunkDecoding end
    struct Unspecified_Chunk <:ChunkDecoding end

    is_DIF_string(s::AbstractString) = occursin(DIF_regexp,s)
    is_DUP_string(s::AbstractString) = occursin(DUP_regexp,s)
    function get_chunk_decoding(s::AbstractString)
        is_DIF = is_DIF_string(s)
        is_DUP = is_DUP_string(s)
        is_DIF && is_DUP && return DIF_DUP
        is_DIF && return DIF
        is_DUP && return DUP
        return No_Chunk_Decoding
    end
    (::Type{DIF})(s::AbstractString) = replace(s,DIF_shift...)
    (::Type{DUP})(s::AbstractString) = replace(s,DUP_shift...)
    (::Type{DIF_DUP})(s::AbstractString) = s |> DIF |> DUP
    (::Type{Unspecified_Chunk})(s::AbstractString) = get_chunk_decoding(s)(s)
    parse_chunk(::Type{No_Chunk_Decoding},s::AbstractString) = Base.parse(Float64,s)
    parse_chunk(::Type{Unspecified_Chunk},s::AbstractString) = parse_chunk(get_chunk_decoding(s),s)
    """
    parse_chunk(::Type{<:ChunkDecoding},s::AbstractString)

Functions to parse data chunk, 
"""
function parse_chunk end
    function parse_chunk(::Type{DUP},s::AbstractString)
        m = match(DUP_regexp,s)
        !isnothing(m) || return Base.parse(Float64,s)
        n = DUP_digits[m.match[1]] # multiplyer numeric
        return  (Base.parse(Float64,s[1:m.offset-1]), n, DUP)
    end
    function parse_chunk(::Type{DIF},s::AbstractString)
        is_DIF_string(s) || return Base.parse(Float64,s)
        val = Base.parse(Float64,replace(s,DIF_digits...))
        return (val,DIF)
    end
    function parse_chunk(::Type{DIF_DUP},s::AbstractString)
        #@show s
        is_DUP = is_DUP_string(s)
        is_DIF = is_DIF_string(s)
        if is_DIF && is_DUP
            m = match(DUP_regexp,s)
            (val,) = parse_chunk(DIF,s[1:m.offset-1])
            n = DUP_digits[m.match[1]] # multiplyer numeric
            return (val,n,DIF_DUP)
        end
        is_DUP && return parse_chunk(DUP,s)
        return parse_chunk(DIF,s)
    end
    function set_or_push!(v::Vector{T}, i::Int, val::T) where T<:Number
        M = length(v)
        if 1 <= i <= M
            v[i] = val
        elseif i == M + 1
            push!(v, val)
        else i > M + 1
            resize!(v, i) 
            v[i] = val
        end
        return 1
    end
    function set_or_push!(v,i,val,::Type{DIF})
        set_or_push!(v,i,v[i-1] + val)
    end
    function set_or_push!(v::AbstractVector{T},i::Int,val,repeates_number,::Type{DIF_DUP}) where {T<:Number}
        N = repeates_number
        N==1 && return set_or_push!(v,i,val) #+ set_or_push!(v,i+1,val)
        M = length(v)
        i+N-1 <= M || resize!(v,i+N-1)
        for jj=0:N-1
            set_or_push!(v,jj+i,val,DIF)
        end
        return N
    end
    function set_or_push!(v::AbstractVector{T},i::Int,val,repeates_number,::Type{DUP}) where {T<:Number}
        N = repeates_number
        N==1 && return set_or_push!(v,i,val,DIF) #+ set_or_push!(v,i+1,val)
        M = length(v)
        i+N-1 <= M || resize!(v,i+N-1)
        vi = @view v[i:i+N-1]
        fill!(vi,val)
        return N
    end
    function set_or_push!(v::MVector{N,T}, i::Int, val::T) where {N,T<:Number}
        v[i]=val
        return 1
    end
"""
    DATAline type specifies the data organization patter viz XY...XY or XY...Y
"""
    abstract type DATAline end
    struct XYYline<:DATAline end
    struct XYXYline<:DATAline end
   
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
        filled_points_counter::Int # stores the number of filled points 
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
    """
    JDXblock(file_name::String)

Creates an empty JDXblock object from full file name 
"""
    function JDXblock(file_name::String) # external constructor
        if !isfile(file_name)
            error("Input filename must be the fullfile name to the existing file")
        end
        jdx = JDXblock()
        jdx.file_name = file_name
        return jdx
    end
"""
    DataBuffer is an itermediate container for data parsed from each string
"""
    mutable struct DataBuffer{  DataLineType <: DATAline, # line format XYYline or XYXYline
                                BufferType   <: AbstractVector,# data buffer type may be static vector or dynamic vector
                                LineType     <: LineDecoding, # this type performs some work on line before parsing 
                                ChunkType    <: ChunkDecoding # this is a type of single chunk to be prased to using Base.parse
                            }
        buffer::BufferType
        xpoints_per_line::IntOrNothing # number of X points in a single line
        ypoints_per_line::IntOrNothing# number of Y points in a single line
        DataBuffer() = new{XYYline,Vector{Float64},No_Line_Decoding,No_Chunk_Decoding}(Float64[],nothing,nothing)
    
        """
        DataBuffer(buffer::BufferType;
                            DataLineType = XYYline,
                            LineDecodingType = No_Line_Decoding,
                            ChunkType = No_Chunk_Decoding,
                            xpoints_per_line::IntOrNothing=nothing,
                            ypoints_per_line::IntOrNothing=nothing) where BufferType<:AbstractVector

    `buffer` - Vector or MVector (if it is known that the number of data hunks is the same for each line)
"""
function DataBuffer(buffer::BufferType;
                            DataLineType = XYYline,
                            LineDecodingType = No_Line_Decoding,
                            ChunkType = No_Chunk_Decoding,
                            xpoints_per_line::IntOrNothing=nothing,
                            ypoints_per_line::IntOrNothing=nothing) where BufferType<:AbstractVector 

            return new{DataLineType,BufferType,LineDecodingType,ChunkType}(buffer,xpoints_per_line,ypoints_per_line)
        end
    end
    is_resizable(db::DataBuffer) = buffertype(db) <:Vector
    has_DIF_type(::DataBuffer{D,B,L,C}) where {D,B,L,C} =  L <: DIF || L <: DIF_DUP 
    datalinetype(::DataBuffer{D}) where D = D
    buffertype(::DataBuffer{D,B}) where {D,B} = B
    linedecodingtype(::DataBuffer{D,B,L}) where {D,B,L} = L
    chunktype(::DataBuffer{D,B,L,C}) where {D,B,L,C} = C
    Base.length(db::DataBuffer{D,B}) where {D,B<:Vector} = length(db.buffer)
    Base.length(::DataBuffer{D,MVector{N,T}}) where {D,N,T} = N

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
    read_jdx_file(file_name::String;
                fixed_columns_number::Bool=false,
                delimiter = nothing,
                validate_data::Bool=true)

Reads JCAMP-DX file `file_name`

Input arguments: 

`file_name` - full file name
(optional keyword args) 
delimiter  - data chunks delimiter (default value is space)
`fixed_columns_number` - if it is known 
that each line in file has the same number of data chunks, this flag 

returns namedtuple (or a vector of namedtuples in the case of multiple blocks) 

with fields:

x - coordinate (wavelength, wavenumber or other)

y - data

headers - dictionary in "String => value" format with headers values 


"""
function read_jdx_file(file_name::String;
                fixed_columns_number::Bool=false,
                delimiter = nothing,
                fixed_line_decoding::Bool = false,
                fixed_chunk_decoding::Bool = false,
                validate_data::Bool=true)

        @assert isfile(file_name) "Must be a file"
        jdx_blocks = find_blocks(file_name)
        if length(jdx_blocks) ==1
            return read!(jdx_blocks[1],
                        delimiter=delimiter,
                        fixed_columns_number=fixed_columns_number,
                        fixed_line_decoding = fixed_line_decoding,
                        fixed_chunk_decoding = fixed_chunk_decoding,
                        validate_data = validate_data)
        else
            return read!.(jdx_blocks,
                            delimiter=delimiter,
                            fixed_columns_number=fixed_columns_number,
                            fixed_line_decoding = fixed_line_decoding,
                            fixed_chunk_decoding = fixed_chunk_decoding,
                            validate_data = validate_data)
        end
    end
    """
    parse_headers!(jdx::JDXblock,headers::Vector{String})

internal function fills headers dictionary from a vector of strings
"""
function parse_headers!(jdx::JDXblock,headers::Vector{String})
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

    """
    addline!( jdx::JDXblock, 
                        data_buffer::DataBuffer{DataLineType,B,LD,ChunkDecoding},
                        current_line::String; # index of current data chunk
                        delimiter=DEFAULT_DELIMITER(DataLineType),
                        validate_data::Bool=true) where {DataLineType<:XYYline,B,LD,ChunkDecoding}

This function parses `current_line` string of file and fills the 
parsed data to the y-vector of `jdx`  object `number_of_y_point_per_chunk` 
    - the number of data point (excluding the x-coordinate) in the line 
    `chunk_index` index of chunk 
    `delimiter` data points delimiter used in `split` function

"""
    function addline!( jdx::JDXblock, 
                        data_buffer::DataBuffer{DataLineType,B,LD,ChunkDecoding},
                        current_line::String; # index of current data chunk
                        delimiter=DEFAULT_DELIMITER(DataLineType),
                        validate_data::Bool=true) where {DataLineType<:XYYline,B,LD,ChunkDecoding}
        
        if isempty(current_line) || current_line[1]=='#'  
             return NaN
        end 
        is_first_line = jdx.filled_points_counter == 0
        # returns the total number of points added to the buffer
        # if the buffer is resizable it is equal to the length of the buffer  
        cur_points_number =  fill_data_buffer!(data_buffer,current_line,delimiter)
        data_buffer.ypoints_per_line = cur_points_number - 1 
        data_buffer.xpoints_per_line = 1

        buffer_starting_index = 2
        buffer_ending_index = cur_points_number
        data_starting_index = 1 + jdx.filled_points_counter
        data_ending_index =  data_starting_index + cur_points_number - 2
        if is_first_line
            resize!(jdx.y_data, cur_points_number - 1 )
        elseif ChunkDecoding <: Union{DIF,DIF_DUP} || 
            (ChunkDecoding <: Unspecified_Chunk &&
             get_chunk_decoding(current_line) <: Union{DIF,DIF_DUP})
            # this works if the chunk decoding type is dif, than the first y value of each line after the first only_headers
            # is used for checks
            buffer_starting_index += 1
            data_ending_index -=1
        end 
        v = @view jdx.y_data[ data_starting_index : data_ending_index ]
        y_data_buffer =@view data_buffer.buffer[buffer_starting_index:buffer_ending_index]
        copyto!(v,y_data_buffer)
        jdx.filled_points_counter += buffer_ending_index - buffer_starting_index + 1
        return data_buffer.buffer[1] # returns x-value for checks
    end
    """
    addline!(jdx::JDXblock, 
                        data_buffer::DataBuffer{DataLineType},
                        current_line::String; # index of current data chunk
                        delimiter=r"[,;]",
                        validate_data::Bool=true) where DataLineType<:XYXYline

Adds line to the XY...XY data
"""
function addline!(jdx::JDXblock, 
                        data_buffer::DataBuffer{DataLineType},
                        current_line::String; # index of current data chunk
                        delimiter=r"[,;]",
                        validate_data::Bool=true) where DataLineType<:XYXYline

        if isempty(current_line) || current_line[1]=='#'  
                    return NaN
        end 
        is_first_line = jdx.filled_points_counter == 0
        cur_points_number =  fill_data_buffer!(data_buffer,current_line,delimiter)
        if is_first_line
            number_of_y_point_per_chunk = div(cur_points_number,2) 
            data_buffer.ypoints_per_line = number_of_y_point_per_chunk
            data_buffer.xpoints_per_line = cur_points_number - number_of_y_point_per_chunk
            resize!(jdx.y_data,number_of_y_point_per_chunk)
            resize!(jdx.x_data,data_buffer.xpoints_per_line)   
        end 
        starting_index = 1 + jdx.filled_points_counter
        ending_index =  starting_index + div(cur_points_number,2) - 1
        vY = @view jdx.y_data[ starting_index : ending_index ]
        vX = @view jdx.x_data[ starting_index : ending_index ]
        y_data_buffer = @view data_buffer.buffer[2:2:cur_points_number]
        x_data_buffer = @view data_buffer.buffer[1:2:cur_points_number]
        copyto!(vY,y_data_buffer)
        copyto!(vX,x_data_buffer)
        jdx.filled_points_counter += div(cur_points_number,2) 
        return data_buffer.buffer[1] # returns x-value for checks
    end   
    """
    fill_data_buffer!(data_buffer::Vector{Float64},current_line,delimiter,chunk_counter::Int=1)

Appends data to vector, initial chunk can be of zero size

"""
    function fill_data_buffer!(data_buffer::DataBuffer{D,B,LineDecodingType,ChunkType},
                                            current_line,
                                            delimiter, 
                                            chunk_counter=1) where {D,B,
                                            LineDecodingType,
                                            ChunkType}
        line_decoded = current_line |> LineDecodingType |> ChunkType
        for chunk in eachsplit(line_decoded,delimiter)
            chunk = strip(chunk)
            !isempty(chunk) || continue # check for empty string
            #chunk_counter += split_data_chunk!(data_buffer,s,chunk_counter)
            chunk_counter +=set_or_push!(data_buffer.buffer,chunk_counter,parse_chunk(ChunkType,chunk)...)
        end
        points_in_chunk = chunk_counter-1
        is_resizable(data_buffer) || length(data_buffer) != points_in_chunk && resize!(data_buffer.buffer,points_in_chunk)
        return points_in_chunk
    end

    """
    generateVectors!(::JDXblock,::Type{<:DATAline})

Generates x-vector and precreates data vectors
"""
function generateVectors!(::JDXblock,::Type{<:DATAline})end
function generateVectors!(jdx::JDXblock,::Type{XYYline})
            point_number = round(Int64,jdx.data_headers["NPOINTS"])
            starting_X = haskey(jdx.data_headers,"FIRSTX") ? jdx.data_headers["FIRSTX"] : 0.0
            step_value =  (jdx.data_headers["LASTX"] -  starting_X)/(point_number-1)
            if haskey(jdx.data_headers,"DELTAX") # check for the delta x values
                check_data_point!(jdx.violation,
                                    jdx.data_headers["DELTAX"],
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

Parses headers from JCAMP-DX file? returns dictionary with file headers
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
                               fixed_columns_number::Bool=true,
                               fixed_line_decoding::Bool = false,
                               fixed_chunk_decoding::Bool = false,
                               validate_data::Bool = true)


fills precreated JDXblock object see [`JDXblock`](@ref).

    `delimiter` - data chunk delimiters
    `only_headers`  - if true parses only block headers
    `fixed_columns_number`  if true, all lines are supposed to have the same number of chunks
    `fixed_line_decoding` if true, data line decoding (No_Line_Decoding, SQZ,PAC) is supposed  to be the same for all lines (parser obtains decoding from the first line of data)
    `fixed_chunk_decoding` if true, all data chunks decoding (No_Chunk_Decoding,DIF,DUP) is supposed  to be the same for all lines (parser obtains decoding from the first line of data)
    `validate_data` if false: we don't need no validation
"""
function read!(jdx::JDXblock; delimiter=nothing,
                               only_headers::Bool=false,
                               fixed_columns_number::Bool=true,
                               fixed_line_decoding::Bool = false,
                               fixed_chunk_decoding::Bool = false,
                               validate_data::Bool = true)
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
                if is_data_started 
                    if !only_headers
                        delimiter =  isnothing(delimiter) ? DEFAULT_DELIMITER(line_type) : delimiter
                        first_line_buffer = DataBuffer(Vector{Float64}(undef,6),
                                                    DataLineType =line_type,
                                                    LineDecodingType = get_line_decoding(ln),
                                                    ChunkType = get_chunk_decoding(ln) )
                                                    
                        x_point = addline!( jdx, first_line_buffer,ln, delimiter = delimiter) # this function modifies
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
            parse_headers!(jdx,header_lines) #parses headers to dict from a vector of strings
            if only_headers # if we need only block headers without data 
                return (x=Vector{Float64}([]),
                        y=Vector{Float64}([]), 
                        headers=jdx.data_headers,
                        is_violated = false)
            end
            number_of_y_point_per_chunk = length(jdx.y_data) # number of numeric points per line
            # by default we assume that the last line is ##END= thus it is ignored
            data_lines_number = total_number_of_lines - header_lines_counter - 1 # last line is ignored because there should be ##END
            if haskey(jdx.data_headers,"NPOINTS") # correct JDX file
                total_point_number = round(Int64,jdx.data_headers["NPOINTS"]) # total number of points is known
                is_XYYline = line_type <: XYYline
                #resize!(jdx.y_data,total_point_number)
                generateVectors!(jdx,line_type)
                #is_XYYline ? generateXvector!(jdx) : resize!(jdx.x_data,total_point_number)
                points_number_per_chunk = length(first_line_buffer)
                # if fix_decoding the line decoding type is taken from the first line, otherwise 
                # line type check is performed each line scan
                line_decoding_type = !fixed_line_decoding ? Unspecified_Line : linedecodingtype(first_line_buffer)
                chunk_decoding_type = !fixed_chunk_decoding ? Unspecified_Chunk : chunktype(first_line_buffer)
                buffer = !fixed_columns_number ? copy(first_line_buffer.buffer) : MVector{points_number_per_chunk,Float64}(undef)
                if !fixed_columns_number && fixed_chunk_decoding && fixed_line_decoding
                    data_buffer = first_line_buffer
                else
                    data_buffer = DataBuffer(buffer,
                                             DataLineType=line_type,
                                             LineDecodingType = line_decoding_type,
                                             ChunkType = chunk_decoding_type,
                                             xpoints_per_line = first_line_buffer.xpoints_per_line,
                                             ypoints_per_line = first_line_buffer.ypoints_per_line)
                end    
                x_factor = get(jdx.data_headers,"XFACTOR",1.0)
                if is_XYYline    
                    x_gen = jdx.x_data[1] # generated x to compare
                    for i in 2:data_lines_number
                        ln = readline(io_file)
                        x_point = addline!(jdx,data_buffer,ln,delimiter=delimiter)
                        x_point = x_point * x_factor
                        if !isnan(x_point) #&& fixed_columns_number
                            x_gen = jdx.x_data[(i-1)*number_of_y_point_per_chunk + 1] # generated x
                        else
                            break
                        end#endif       
                    end#endfor
                else #XYXYline
                    for i in 2:data_lines_number
                        ln = readline(io_file)
                        addline!(jdx,data_buffer,ln,delimiter=delimiter)     
                    end
                end
            else number_of_y_point_per_chunk==1 # file with two columns like CSV
                total_point_number = data_lines_number
                resize!(jdx.y_data,total_point_number)
                jdx.x_data = similar(jdx.y_data)
                jdx.x_data[1] = x_point
                data_buffer = DataBuffer(MVector{2,Float64}(undef))
                for i in 2:data_lines_number
                    ln = readline(io_file)
                     x_out = addline!(jdx,data_buffer, ln,delimiter=delimiter)
                     !isnan(x_out) ? jdx.x_data[i] = x_out : nothing
                end
            end

        end # close file
        if haskey(jdx.data_headers,"YFACTOR") 
            y_factor::Float64 = jdx.data_headers["YFACTOR"] 
            jdx.y_data .*= y_factor
        end        
        return (x=jdx.x_data,
                y=jdx.y_data, 
                headers=jdx.data_headers)
    end#endof read!

    """
    find_blocks(file_name)

Counts blocks in file `file_name`, fills coordinates of blocks start
and finish in file lines, returns the vector of [`JDXblock`](@ref) blocks
or a singe block.
"""
function find_blocks(file_name)
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
end
