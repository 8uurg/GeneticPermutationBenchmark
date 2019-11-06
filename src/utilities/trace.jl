
## Some progress tracing helpers

struct ProgressTrace
    # Starting time.
    starting_time :: Ref{Float64}
    # Best value found so far.
    best :: Ref{Float64}
    # Vector to store in results
    results :: Vector{Float64}
    # Time up to which the result should be accounted for (for a certain bin)
    moments :: Vector{Float64}
    cmoment :: Ref{Int64}
    # Target
    target :: Float64
    # Hitting time
    hitting_time :: Ref{Float64}

    function ProgressTrace(moments :: Vector{Float64}, target :: Float64)
        new(time(), Ref(0.0), zeros(length(moments)), sort(moments), Ref(1), target, Ref(typemax(Float64)))
    end
end

function update_trace!(trace :: ProgressTrace, objective :: Float64)
    if trace.best[] >= objective
        return nothing
    end
    trace.best[] = objective
    tdiff = time() - trace.starting_time[]
    #println("Improvement at: $tdiff (s), new value is: $(objective)")
    while trace.cmoment[] <= length(trace.results) && trace.moments[trace.cmoment[]] < tdiff
        trace.cmoment[] += 1
        if trace.cmoment[] <= length(trace.results)
            trace.results[trace.cmoment[]] = trace.results[trace.cmoment[] - 1]
        end
    end
    if trace.cmoment[] <= length(trace.results)
        trace.results[trace.cmoment[]] = objective
        #println("tdiff ($tdiff) < $(trace.moments[trace.cmoment[]]), updating results[$(trace.cmoment[])] = $(objective).")
    end
    if objective >= trace.target && trace.hitting_time[] > tdiff
        trace.hitting_time[] = tdiff
        #println("Hit Optimum at $(tdiff) with $(objective)!")
    end
    return nothing
end

function update_trace!(:: Nothing, objective :: Float64)
end

function postprocess_trace!(trace :: ProgressTrace)
    cb = trace.results[1]
    for i in 1:length(trace.results)
        cb = max(cb, trace.results[i])
        trace.results[i] = cb
    end
end

function postprocess_trace!(::Nothing)
end

function trace_bb(f :: Function, trace :: Union{Nothing,ProgressTrace})
    function rex(assignment :: Vector{Int64})
        fitness = f(assignment)
        update_trace!(trace, convert(Float64, fitness))
        return fitness
    end
    return rex
end