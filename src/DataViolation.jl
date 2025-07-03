 """
        Stores violation points
    """
    struct ViolationPoint{T<:Number}
        value_checker::T
        value_real::T
        difference::T
        point_index::Int
        line_index::Int
        message::String
    end
    const X_VIOLATION_CRITERIUM = 1e-2
    const DELTAX_VIOLATION_CRITERIUM = 1e-3
    const Y_VIOLATION_CRITERIUM = 1e-2
mutable struct DataViolation
        x::Union{Vector{ViolationPoint},Nothing}
        y::Union{Vector{ViolationPoint},Nothing}
        npoints::Union{ViolationPoint,Nothing}
        deltax::Union{ViolationPoint,Nothing}
        DataViolation() = new(nothing,nothing,nothing,nothing)
end

function Base.push!(dv::DataViolation,p::ViolationPoint,field::Symbol)
    is_vector = field==:x || field ==:y
    if !is_vector
        setfield!(dv,field,p)
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

function check_data_point!(dv::DataViolation,value_checker::T,
                                    value_real::T,field::Symbol,
                                    point_index::Int=0,
                                    line_index::Int=0) where T
    if field == :x
        criterium = X_VIOLATION_CRITERIUM
    elseif field == :y
        criterium = Y_VIOLATION_CRITERIUM
    elseif field == :deltax
        criterium = DELTAX_VIOLATION_CRITERIUM
    elseif field == :npoints
        criterium = 1.0
    else
        error("Unknown violation checker fieldname")
    end
    difference = abs(value_checker - value_real)
    if  !isapprox(value_real, value_checker,rtol=criterium)
        message = "$(field) is violated $(point_index>0 ? "at point"*string(point_index) : nothing), 
                    $(line_index>0 ? "at line"*string(line_index) : nothing), 
                    the value of approx(val_calculated,val_checker,rtol= $(criterium))"
        push!(dv,ViolationPoint(value_checker,value_real,difference,point_index,line_index,message))
    end

end