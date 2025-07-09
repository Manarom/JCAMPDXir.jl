 """
        Stores Validation points
    """
    struct ValidationPoint{T<:Number}
        value_checker::T
        value_real::T
        difference::T
        point_index::Int
        line_index::Int
        message::String
    end
    Base.length(::ValidationPoint) = 1
    const ValpNothing = Union{ValidationPoint,Nothing}
    const X_VIOLATION_CRITERIUM = 1e-2
    const DELTAX_VIOLATION_CRITERIUM = 1e-3
    const Y_VIOLATION_CRITERIUM = 1e-2

    const VIOLATION_CRITERIA = Dict( :npoints => 1.0,
                                     :deltax => DELTAX_VIOLATION_CRITERIUM ,
                                     :maxx => X_VIOLATION_CRITERIUM ,
                                     :minx => X_VIOLATION_CRITERIUM ,
                                     :miny => Y_VIOLATION_CRITERIUM ,
                                     :maxy => Y_VIOLATION_CRITERIUM )
    mutable struct DataValidation
        x::Union{Vector{ValidationPoint},Nothing}
        y::Union{Vector{ValidationPoint},Nothing}
        single_point_validators::Dict{Symbol,ValpNothing}
        DataValidation() = begin
            new(nothing,nothing,Dict{Symbol,ValpNothing}(k=>nothing for k in keys(VIOLATION_CRITERIA)))
        end 
end
Base.show(dv::DataType) = print(validation_message(dv))
function validation_message(dv::DataValidation)
    warning_message = ""
    has_violations(dv) || return warning_message
    if !isnothing(dv.x) 
        warning_message = warning_message*" data check failed for x at $(length(dv.x)) points \n"
    end
    if    !isnothing(dv.y) 
        warning_message = warning_message*" data check failed for y at $(length(dv.y)) points \n"
    end
    for (k,v) in dv.single_point_validators
        !isnothing(v) || continue
        warning_message *= " data check failed for $(k) \n"
    end
    return warning_message
end
function has_violations(dv::DataValidation)
    isnothing(dv.x) && isnothing(dv.y) ? nothing : return true
    for p in values(dv.single_point_validators)
        isnothing(p) || return true
    end
    return false
end
function Base.push!(dv::DataValidation,p::ValidationPoint,field::Symbol)
    is_vector = field==:x || field ==:y
    if !is_vector && haskey(dv.single_point_validators,field)
        setindex!(dv.single_point_validators,p,field)
        return dv
    end
    data = getfield(dv,field)
    if isnothing(data)
        data = fill(p,1)
    else
        push!(data,p)
    end
    return dv
end

function check_data_point!(dv::DataValidation,
                                    value_checker::T,
                                    value_real::T,
                                    field::Symbol,
                                    point_index::Int=0,
                                    line_index::Int=0) where T
    if field == :x
        criterium = X_VIOLATION_CRITERIUM
    elseif field == :y
        criterium = Y_VIOLATION_CRITERIUM
    elseif haskey(dv.single_point_validators,field)
        criterium = VIOLATION_CRITERIA[field]
    else
        error("Unknown Validation checker fieldname")
    end
    difference = abs(value_checker - value_real)
    if  !isapprox(value_real, value_checker,rtol=criterium)
        message = "$(field) is violated $(point_index>0 ? "at point"*string(point_index) : ""), $(line_index>0 ? "at line"*string(line_index) : ""), 
                    the approx($(field)_checker, $(field)_checker) = approx($(value_real),$(value_checker),rtol= $(criterium)) returns false"
        push!(dv,ValidationPoint(value_checker,value_real,difference,point_index,line_index,message),field)
        return false
    end
    return true
end