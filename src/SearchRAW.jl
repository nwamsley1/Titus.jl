function SearchRAW(
                   spectra::Arrow.Table, 
                   ptable::PrecursorDatabase,
                   selectTransitions, 
                   right_precursor_tolerance::Float32,
                   left_precursor_tolerance::Float32,
                   transition_charges::Vector{UInt8},
                   transition_isotopes::Vector{UInt8},
                   b_start::Int64,
                   y_start::Int64,
                   fragment_match_ppm::Float32,
                   ms_file_idx::UInt32,
                   )
    
    scored_PSMs = makePSMsDict(FastXTandem())
    #not goinig to want this for all methods
    allFragmentMatches = Vector{FragmentMatch}()

    #precursorList needs to be sorted by precursor MZ. 
    #Iterate through rows (spectra) of the .raw data. 
    i = 0
    skip = 1
    for spectrum in Tables.namedtupleiterator(spectra)
        i += 1
        if isequal(spectrum[:precursorMZ], missing)
            skip += 1
            continue
        end
        #
        #params = getSpectrumSpecificParams(spectrum, selectParams)

        transitions = selectTransitions(spectrum[:precursorMZ], 
                                        ptable,
                                        right_precursor_tolerance,
                                        left_precursor_tolerance
                                        )
        #Named tuple for scan 
        fragmentMatches = matchPeaks(transitions, spectrum[:masses], spectrum[:intensities], 
                                    #δs = params[:δs],
                                    δs = Float64[0.0],#[Float64(0)],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx)
        #sort(getTransitions(selectPrecursors(spectrum)),  
        #                    by = t-> getMZ(t))
        unscored_PSMs = UnorderedDictionary{UInt32, FastXTandem}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)
        
        #Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(spectrum[:scanNumber]))
        Score!(scored_PSMs, unscored_PSMs, scan_idx = Int64(i))

        append!(allFragmentMatches, fragmentMatches)
    end
    #println(prec_id_to_transitions)
    (scored_PSMs, allFragmentMatches)
end