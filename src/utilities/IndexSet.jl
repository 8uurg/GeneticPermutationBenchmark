## Index Set related utilities
# This nicely outperforms the standard BitSet in this context.

using StaticArrays
# import Base:push!,pop!,append!,iterate,length,setdiff,issubset,isless,==,hash,&

# function _div64(v :: Int64)
#     div(v, 64) #v >> 6 should do the same thing
# end
# function _mod64(v :: Int64)
#     mod(v, 64)
# end
import Base._div64
import Base._mod64

function _setbit(nt :: UInt64, v :: Int64)
    nt | UInt64(1) << (v)
end
function _setbit(nt :: UInt64, v :: UInt64)
    nt | UInt64(1) << (v)
end

function _unsetbit(nt :: UInt64, v :: Int64)
    nt & ~(UInt64(1) << (v))
end
function _unsetbit(nt :: UInt64, v :: UInt64)
    nt & ~(UInt64(1) << (v))
end

struct IndexSet{S}
    bits :: MVector{S, UInt64}

    function IndexSet(n :: Int64)
        s = _div64(n)+1
        new{s}(MVector{s}(zeros(UInt64, s)))
    end

    function IndexSet(o :: IndexSet{S}) where {S}
        new{S}(copy(o.bits))
    end
end

function Base.copy(s :: IndexSet{S}) :: IndexSet{S} where {S}
    return IndexSet(s)
end

function Base.length(os :: IndexSet)
    sum(count_ones.(os.bits))
end

function Base.push!(os :: IndexSet, v :: UInt64)
    @inbounds os.bits[_div64(v)+1] = _setbit(os.bits[_div64(v)+1], _mod64(v))
    os
end

function Base.push!(os :: IndexSet, v :: Int64)
    @inbounds os.bits[_div64(v)+1] = _setbit(os.bits[_div64(v)+1], _mod64(v))
    os
end

function Base.push!(os :: IndexSet{S}, v :: IndexSet{S}) where {S}
    os.bits .|= v.bits
    os
end

function Base.pop!(os :: IndexSet, v :: Int64)
    @inbounds os.bits[_div64(v)+1] = _unsetbit(os.bits[_div64(v)+1], _mod64(v))
    os
end

function Base.empty!(set :: IndexSet)
    set.bits .= UInt64(0)
end

function Base.copyto!(dest :: IndexSet, from :: IndexSet)
    copyto!(dest.bits, from.bits)
end

function Base.append!(os :: IndexSet, vs :: IndexSet)
    os.bits .|= vs.bits
    os
end

function Base.:&(as :: IndexSet, bs :: IndexSet)
    l = min(length(as.bits), length(bs.bits))
    bits = zeros(UInt, l)
    for i in 1:l
        @inbounds bits[i] = as.bits[i] & bs.bits[i]
    end
    IndexSet(bits)
end

# Fast-hash
function Base.hash(os :: IndexSet, h :: UInt=zero(UInt))
    if length(os.bits) == 0
       return hash(h)
    end
    for bit in os.bits
        h ⊻= hash(bit)
    end
    return hash(h) # Maybe call hash here as well.
end

# Note Base.:== does not work anymore for some reason.
# Workaround.
import Base.==
function ==(a :: IndexSet, b :: IndexSet)
    if length(a.bits) > length(b.bits)
        for abit in view(a.bits, length(b.bits)+1, length(a.bits))
            if !iszero(abit)
                return false
            end
        end
    elseif length(b.bits) > length(a.bits)
        for bbit in view(b.bits, length(a.bits)+1, length(b.bits))
            if !iszero(bbit)
                return false
            end
        end
    end
    for (abit, bbit) in zip(a.bits, b.bits)
        if abit != bbit
            return false
        end
    end
    return true
end

function _togglebit!(os :: IndexSet, i :: Int64)
    msk = UInt64(1) << _mod64(i)
    @inbounds os.bits[div(i, 64)+1] ⊻= msk
    return @inbounds (os.bits[div(i, 64)+1] & msk) == msk
end

function _getbit(os :: IndexSet, i :: Int64) :: Bool
    msk = UInt64(1) << _mod64(i)
    @inbounds (os.bits[div(i, 64)+1] & msk) == msk
end

function Base.iterate(os :: IndexSet{S}) :: Union{Nothing, Tuple{Int64, Tuple{Int64, UInt64}}} where {S}
    # Empty set, no memory initialized.
    if length(os.bits) == 0
        return nothing
    end
    statei = 1
    statebits = os.bits[1]
    iterate(os, (statei, statebits))
end

function Base.iterate(os :: IndexSet{S}, state :: Tuple{Int64, UInt64}) :: Union{Nothing, Tuple{Int64, Tuple{Int64, UInt64}}} where {S}
    statei, statebits = state
    tz = trailing_zeros(statebits)
    while tz == 64 && statei < length(os.bits)
        statei += 1
        @inbounds statebits = os.bits[statei]
        tz = trailing_zeros(statebits)
    end
    if tz == 64
        return nothing
    end
    statebits = _unsetbit(statebits, convert(UInt64, tz))
    (tz + (statei - 1) * 64, (statei, statebits))
end

# Lazy set difference for iteration.
struct SetDiff{S}
    a :: IndexSet{S}
    b :: IndexSet{S}
end

function Base.setdiff!(a :: IndexSet{S}, b :: IndexSet{S}) where {S}
    a.bits .&= .~ b.bits
    a
end

function lazysetdiff(a :: IndexSet{S}, b :: IndexSet{S}) where {S}
    SetDiff{S}(a, b)
end

function Base.iterate(sd :: SetDiff{S}) where {S}
    if length(sd.a.bits) == 0 || length(sd.b.bits) == 0
        return nothing
    end
    statei = 1
    @inbounds statebits = sd.a.bits[1] & ~sd.b.bits[1]
    iterate(sd, (statei, statebits))
end

function Base.iterate(sd :: SetDiff{S}, state :: Tuple{Int64, UInt64}) :: Union{Nothing, Tuple{Int64, Tuple{Int64, UInt64}}} where {S}
    statei, statebits = state
    tz = trailing_zeros(statebits)
    while tz == 64 && statei < min(length(sd.a.bits), length(sd.b.bits))
        statei += 1
        @inbounds statebits = sd.a.bits[statei] & ~sd.b.bits[statei]
        tz = trailing_zeros(statebits)
    end
    if tz == 64
        return nothing
    end
    statebits = _unsetbit(statebits, tz)
    (tz + (statei - 1) * 64, (statei, statebits))
end

# Lazy intersection for IndexSets
struct SetIntersection{S}
    a :: IndexSet{S}
    b :: IndexSet{S}
end

function Base.intersect!(a :: IndexSet{S}, b :: IndexSet{S}) where {S}
    a.bits .&= b.bits
    a
end

function lazyintersect(a :: IndexSet{S}, b :: IndexSet{S}) where {S}
    return SetIntersection(a, b)
end

function Base.iterate(sd :: SetIntersection{S}) where {S}
    if length(sd.a.bits) == 0 || length(sd.b.bits) == 0
        return nothing
    end
    statei = 1
    @inbounds statebits = sd.a.bits[1] & sd.b.bits[1]
    iterate(sd, (statei, statebits))
end

function Base.iterate(sd :: SetIntersection, state :: Tuple{Int64, UInt64}) :: Union{Nothing, Tuple{Int64, Tuple{Int64, UInt64}}}
    statei, statebits = state
    tz = trailing_zeros(statebits)
    while tz == 64 && statei < min(length(sd.a.bits), length(sd.b.bits))
        statei += 1
        @inbounds statebits = sd.a.bits[statei] & sd.b.bits[statei]
        tz = trailing_zeros(statebits)
    end
    if tz == 64
        return nothing
    end
    statebits = _unsetbit(statebits, tz)
    (tz + (statei - 1) * 64, (statei, statebits))
end

function Base.issubset(a::IndexSet, b::IndexSet)
    # `a` is not a subset if there are additional bits that are not zero.
    if length(a.bits) > length(b.bits)
        for i in length(b.bits)+1:length(a.bits)
            if !iszero(a.bits[i])
                return false
            end
        end
    end
    # Otherwise, masking a with b, should yield a if a is a subset.
    for (ab, bb) in zip(a.bits, b.bits)
        if ab & bb != ab
            return false
        end
    end
    return true
end
