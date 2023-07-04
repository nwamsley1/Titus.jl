function findFirstFragmentBin(frag_index::Vector{FragBin{T}}, frag_min::AbstractFloat, frag_max::AbstractFloat) where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    potential_match = nothing
    while lo <= hi

        mid = (lo + hi) ÷ 2

        if (frag_min) <= getHighMZ(frag_index[mid])
            if (frag_max) >= getHighMZ(frag_index[mid]) #Frag tolerance overlaps the upper boundary of the frag bin
                potential_match = mid
            end
            hi = mid - 1
        elseif (frag_max) >= getLowMZ(frag_index[mid]) #Frag tolerance overlaps the lower boundary of the frag bin
            if (frag_min) <= getLowMZ(frag_index[mid])
                potential_match = mid
                #return mid
            end
            lo = mid + 1
        end
    end

    return potential_match#, Int64(getPrecBinID(frag_index[potential_match]))
end

function searchPrecursorBin!(precs::Counter{UInt32, UInt8, Float32}, intensity::Float32, precursor_bin::PrecursorBin{T}, window_min::Float64, window_max::Float64) where {T<:AbstractFloat}
   
    N = getLength(precursor_bin)
    lo, hi = 1, N

    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) < window_min
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    window_start = (lo <= N ? lo : return)

    if getPrecMZ(getPrecursor(precursor_bin, window_start)) > window_max
        return 
    end

    lo, hi = window_start, N

    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getPrecMZ(getPrecursor(precursor_bin, mid)) > window_max
            hi = mid - 1
        else
            lo = mid + 1
        end
    end

    window_stop = hi

    function addFragmentMatches!(precs::Counter{UInt32, UInt8, Float32}, intensity::Float32, precursor_bin::PrecursorBin{T}, start::Int, stop::Int) where {T<:AbstractFloat}
        #@inbounds @simd for precursor_idx in start:stop
        for precursor_idx in start:stop
            prec = getPrecursor(precursor_bin, precursor_idx)
            inc!(precs, getPrecID(prec), getIntensity(prec), intensity)
        end

        return
    end

    addFragmentMatches!(precs, intensity, precursor_bin, window_start, window_stop)

    return 

end

function queryFragment!(precs::Counter{UInt32, UInt8, Float32}, intensity::Float32, frag_index::FragmentIndex{T}, min_frag_bin::Int64, frag_min::Float64, frag_max::Float64, prec_mz::Float32, prec_tol::Float64) where {T<:AbstractFloat}
    
    frag_bin = findFirstFragmentBin(getFragBins(frag_index), frag_min, frag_max)
    #No fragment bins contain the fragment m/z
    if (frag_bin === nothing)
        return min_frag_bin
    #This frag bin has already been searched
    elseif frag_bin <= min_frag_bin
        return min_frag_bin
    end

    i = 1
    prec_min = prec_mz - prec_tol
    prec_max = prec_mz + prec_tol
    #=println("prec_min $prec_min")
    println("prec_max  $prec_max")
    println("precs.size ", precs.size)
    println("frag_min", frag_min)
    println("frag_max ", frag_max)=#
    while (frag_bin < length(getFragBins(frag_index))) #getLowMZ(getFragmentBin(frag_index, frag_bin)) <frag_max
    
        #Fragment bin matches the fragment ion
        #println(i)
        i += 1
        if (getLowMZ(getFragmentBin(frag_index, frag_bin)) > frag_max)
            return frag_bin
        else
            #println("frag_bin $frag_bin")
            #prec_min = prec_mz - prec_tol - 3.0*NEUTRON #These lines are a 25% increase in time. 
            #prec_max = prec_mz + prec_tol + 1.0*NEUTRON
            searchPrecursorBin!(precs, intensity, getPrecursorBin(frag_index, UInt32(frag_bin)), prec_min, prec_max)
            frag_bin += 1
        end

    end

    #Only reach this point if frag_bin exceeds length(frag_index)
    return frag_bin - 1
end

function searchScan!(precs::Counter{UInt32, UInt8, Float32}, f_index::FragmentIndex{T}, min_intensity::U, masses::Vector{Union{Missing, U}}, intensities::Vector{Union{Missing, U}}, precursor_window::AbstractFloat, ppm::AbstractFloat, width::AbstractFloat; topN::Int = 20, min_frag_count::Int = 3) where {T,U<:AbstractFloat}
    
    getFragTol(mass::U, ppm::AbstractFloat) = mass*(1 - ppm/1e6), mass*(1 + ppm/1e6)

    function filterPrecursorMatches!(precs::Counter{UInt32, UInt8, Float32}, topN::Int, min_frag_count::Int) where {T<:AbstractFloat}
        match_count = countFragMatches(precs, min_frag_count)
        prec_count = getSize(precs) - 1
        sort!(precs, topN);
        #precs.size = max(exceeds_min_frag_count, 1)
        return match_count, prec_count#[first(x) for x in sort(filter(x->last(x)>=min_frag_count,collect(precs)), by=x->last(x), rev = true)][1:min(topN, end)], match_count, prec_count
    end
    #println("TEST")
    min_frag_bin = 0
    #println("start")
    #@time begin
    #    println("start ")
        #for i in peak_iterator
        for (mass, intensity) in zip(masses, intensities)

            mass, intensity = coalesce(mass, 0.0),  coalesce(intensity, 0.0)

            if intensity<min_intensity
                continue
            end

            FRAGMIN, FRAGMAX = getFragTol(mass, ppm) 

            min_frag_bin = queryFragment!(precs, intensity, f_index, min_frag_bin, FRAGMIN, FRAGMAX, precursor_window, width)
        end 
    #end
    #println("stop ", length(precs))
    #println("PRECS $precs")
    #println("stop")
    return filterPrecursorMatches!(precs, topN, min_frag_count)
end
