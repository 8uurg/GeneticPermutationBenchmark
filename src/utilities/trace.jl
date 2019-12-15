
## Some progress tracing helpers

struct ProgressTrace
    # Starting time.
    starting_time :: Ref{Float64}
    # Best value found so far.
    best :: Ref{Float64}
    n_evals :: Ref{Int64}
    # Vector to store in results
    results :: Vector{Float64}
    results_eval :: Vector{Float64}
    # Time up to which the result should be accounted for (for a certain bin)
    moments :: Vector{Float64}
    cmoment :: Ref{Int64}
    # Evaluations & index of current bin.
    moments_n_eval :: Vector{Int64}
    cmoment_n_eval :: Ref{Int64}
    # Target
    target :: Float64
    # Hitting time
    hitting_time :: Ref{Float64}
    hitting_eval :: Ref{Int64}

    function ProgressTrace(moments :: Vector{Float64}, target :: Float64)
        new(# Starting Time
            time(), 
            # Best so far
            Ref(typemin(Float64)), 
            Ref(0),
            # Results vectors
            fill(typemin(Float64), length(moments)), # Time
            fill(typemin(Float64), length(moments)), # # evaluations
            # Time moments to measure & current moment
            sort(moments), Ref(1), 
            # Evaluation points to measure & current, empty for compat.
            [], Ref(1),
            # Target
            target, 
            # Hitting time.
            Ref(typemax(Float64)),
            Ref(typemax(Float64)))
    end

    function ProgressTrace(moments :: Vector{Float64}, moments_n_eval :: Vector{Int64}, target :: Float64)
        @assert issorted(moments)
        @assert issorted(moments_n_eval)
        new(# Starting Time
            Ref(time()), 
            # Best so far
            Ref(typemin(Float64)), 
            Ref(0),
            # Results vectors
            fill(typemin(Float64), length(moments)), # Time
            fill(typemin(Float64), length(moments)), # # evaluations
            # Time moments to measure & current moment
            copy(moments), Ref(1), 
            # Evaluation points to measure & current, empty for compat.
            copy(moments_n_eval), Ref(1),
            # Target
            target, 
            # Hitting time.
            Ref(typemax(Float64)),
            Ref(typemax(Int64)))
    end
end

function update_trace!(trace :: ProgressTrace, objective :: Float64)
    # Another evaluation happened.
    trace.n_evals[] += 1
    # Quickly exit on non-improving evaluations.
    if trace.best[] >= objective
        return nothing
    end
    # Solution improved, update trace statistics.
    trace.best[] = objective
    tdiff = time() - trace.starting_time[]
    #println("Improvement at: $tdiff (s), new value is: $(objective)")
    # Seek to next bin for time measurement.
    while trace.cmoment[] <= length(trace.results) && trace.moments[trace.cmoment[]] < tdiff
        trace.cmoment[] += 1
        if trace.cmoment[] <= length(trace.results)
            trace.results[trace.cmoment[]] = trace.results[trace.cmoment[] - 1]
        end
    end
    # Update found bin. If it still exists that is.
    if trace.cmoment[] <= length(trace.results)
        trace.results[trace.cmoment[]] = objective
        #println("tdiff ($tdiff) < $(trace.moments[trace.cmoment[]]), updating results[$(trace.cmoment[])] = $(objective).")
    end
    # Do the same for number of evaluations. Find bin
    while trace.cmoment_n_eval[] <= length(trace.results_eval) && trace.moments_n_eval[trace.cmoment_n_eval[]] < tdiff
        trace.cmoment_n_eval[] += 1
        if trace.cmoment_n_eval[] <= length(trace.results_eval)
            trace.results_eval[trace.cmoment_n_eval[]] = trace.results_eval[trace.cmoment_n_eval[] - 1]
        end
    end
    # Update found bin. If it still exists that is.
    if trace.cmoment_n_eval[] <= length(trace.results_eval)
        trace.results_eval[trace.cmoment_n_eval[]] = objective
    end
    # Update hitting time statistics.
    if objective >= trace.target && trace.hitting_time[] > tdiff
        trace.hitting_time[] = tdiff
        trace.hitting_eval[] = trace.n_evals[]
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