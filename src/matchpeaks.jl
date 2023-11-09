#This function name needs to be more generic?
function rangeQuerySorted(sorted_array::Vector{T}, l_bnd::T, u_bnd::T) where T <: Number
    start = searchsortedfirst(sorted_array, l_bnd ,lt=(t,x)->t<x)
    stop = searchsortedlast(sorted_array, u_bnd,lt=(x, t)->t>x)
    return start, stop
end

function getMS1Peaks!(ms1_max_heights::UnorderedDictionary{UInt32, V}, 
                        precursors::Dictionary{UInt32, Precursor{T}}, 
                        MS1::AbstractArray{Union{Missing, P}}, 
                        INTENSITIES::AbstractArray{Union{Missing, P}}, 
                        precursor_rts::Vector{R}, 
                        precursor_idxs::Vector{UInt32}, 
                        precursor_ms_file_idxs::Vector{UInt32}, 
                        rt::R, rt_tol::R, left_mz_tol::U, right_mz_tol::U, ms_file_idx::UInt32) where {T,P,R,U,V <: AbstractFloat}
    
    #Get precursors for which the best scan RT is within `rt_tol` of the current scan `rt`
    start::Int, stop::Int = rangeQuerySorted(precursor_rts,rt - rt_tol, rt + rt_tol)
    if (stop-start) >= 0
    #Check to see if the MS1 height for each precursor is greater than the maximum previously observed. 
        for (i, psm_idx) in enumerate(start:stop)

            if precursor_ms_file_idxs[psm_idx]!=ms_file_idx
                continue
            end

            precursor_idx = precursor_idxs[psm_idx]
            
            if !isassigned(ms1_max_heights, precursor_idx) #If this precursor has not been encountered before. 
                insert!(ms1_max_heights, precursor_idx, Float32(0))
            end

            mz = getMZ(precursors[precursor_idx]) #Get the precursor MZ

            idx = binaryGetNearest(MS1, mz, mz-left_mz_tol, mz+right_mz_tol) #Get the peak index of the peak nearest in mz to the precursor. 

            if idx == 0
                continue
            end

            if INTENSITIES[idx]>=ms1_max_heights[precursor_idx] #Replace maximum observed MS1 height for the precursor if appropriate. 
                ms1_max_heights[precursor_idx] = INTENSITIES[idx]
            end
        end
    end
end

function getMS1PeakHeights!(ms1_max_heights::UnorderedDictionary{UInt32, T},
                            ptable::PrecursorDatabase, 
                            retentionTimes::AbstractArray, 
                            masses::AbstractArray,
                            intensities::AbstractArray, 
                            msOrders::AbstractArray,
                            precursor_rts::Vector{R}, precursor_idxs::Vector{UInt32}, precursor_ms_file_idxs::Vector{UInt32},
                            rt_tol::R, left_mz_tol::U, right_mz_tol::U, ms_file_idx::UInt32) where {T,R,U<:AbstractFloat}
    #println("tunction")
    #i = 1

    for scan_idx in eachindex(retentionTimes)
        if msOrders[scan_idx]!=Int32(1) #Skip non MS1 scans. 
            continue
        end
        #i += 1
        #println("test?")
        @inline getMS1Peaks!(ms1_max_heights,
                    getIDToPrec(ptable), 
                    masses[scan_idx], 
                    intensities[scan_idx],  
                    precursor_rts, 
                    precursor_idxs, 
                    precursor_ms_file_idxs,                    #Precursor retention times and id's (sorted by rt)                                      #PrecursorTable
                    retentionTimes[scan_idx],                               #RT of current scan
                    rt_tol, left_mz_tol, right_mz_tol, ms_file_idx #Only consider precursors where the best scan is within the `rt_col` of the current scan
                    );
    end
    #println("i ", i)
end