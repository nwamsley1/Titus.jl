function SearchRAW(
                    spectra::Arrow.Table, 
                    #ptable::PrecursorDatabase,
                    frag_index::FragmentIndex{Float32},
                    fragment_list::Vector{Vector{LibraryFragment{Float32}}},
                    ms_file_idx::UInt32,
                    interpolation::Any,
                    err_dist::Laplace{Float64},;
                    isolation_width::Float64 = 4.25,
                    precursor_tolerance::Float64 = 5.0,
                    fragment_tolerance::Float64 = 20.0,
                    frag_ppm_err::Float64 = 0.0,
                    topN::Int64 = 20,
                    min_frag_count::Int64 = 4,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    λ::Float32 = Float32(1e3),
                    γ::Float32 = zero(Float32),
                    scan_range::Tuple{Int64, Int64} = (0, 0), 
                    max_peaks::Int = 200, 
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0),
                    rt_tol::Float32 = Float32(30.0),
                    #fragment_match_ppm::U,
                    data_type::Type{T} = Float64,
                    sample_rate::Float64 = 1.0,
                    collect_frag_errs = false,
                    spec_order::Set{Int64} = Set(2)
                    ) where {T<:Real}
    
    scored_PSMs = makePSMsDict(XTandem(data_type))

    msms_counts = Dict{Int64, Int64}()
    precs = Counter(UInt32, UInt8, Float32, length(fragment_list)) #Prec counter
    n = 1
    all_fmatches = Vector{FragmentMatch{Float32}}()
    collect_matches ? all_fmatches = [FragmentMatch{Float32}() for x in range(1, 30000)] : nothing
    #for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    fragmentMatches = [FragmentMatch{Float32}() for _ in range(1, 100000)]
    fragmentMisses = [FragmentMatch{Float32}() for _ in range(1, 100000)]
    transitions = [LibraryFragment{Float32}() for _ in range(1, 100000)]
    prec_ids = [zero(UInt32) for _ in range(1, 10000)]

    H_COLS = zeros(Int64, 100000)
    H_ROWS = zeros(Int64, 100000)
    H_VALS = zeros(Float32, 100000)
  
    prec_idx = 0
    transition_idx = 0

    minimum_rt, maximum_rt = minimum(interpolation), maximum(interpolation)
    #for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    for i in range(1, size(spectra[:masses])[1])

        (i%10000) == 0 ? println(i) : nothing
        msn = spectra[:msOrder[i]]
        msn ∈ spec_order ? nothing : continue
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1 

        i ∈ scan_range ? nothing : continue
        first(rand(1)) <= sample_rate ? nothing : continue

        min_intensity = zero(Float32)# spectrum[:intensities][sortperm(spectrum[:intensities], rev = true)[min(max_peaks, length(spectrum[:intensities]))]]


        iRT = interpolation(spectra[:retentionTime][i])::Float64
        (iRT < minimum_rt) ? iRT_low = Float32(-Inf) : iRT_low = Float32(iRT - rt_tol)
        (iRT > maximum_rt) ? iRT_high = Float32(Inf) : iRT_high = Float32(iRT + rt_tol)

        prec_count, match_count = searchScan!(precs,
                    frag_index, 
                    min_intensity, spectra[:masses][i], spectra[:intensities][i], spectra[:precursorMZ][i], 
                    iRT_low, iRT_high,
                    Float32(fragment_tolerance), 
                    Float32(precursor_tolerance),
                    Float32(isolation_width),
                    min_frag_count = min_frag_count, 
                    min_ratio = Float32(min_matched_ratio),
                    topN = topN
                    )
        
        transition_idx, prec_idx = selectTransitions!(transitions, fragment_list, precs, topN)
        transition_idx < 2 ? continue : nothing 
        
        nmatches, nmisses = matchPeaks(transitions, 
                                    transition_idx,
                                    fragmentMatches,
                                    fragmentMisses,
                                    spectra[:masses][i], 
                                    spectra[:intensities][i], 
                                    FragmentMatch{Float32},
                                    count_unmatched=true,
                                    δs = [frag_ppm_err],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity,
                                    ppm = fragment_tolerance
                                    )

        if nmatches < 2
            reset!(fragmentMatches, nmatches), reset!(fragmentMisses, nmisses)
            continue
        end

        n = collectFragErrs(all_fmatches, fragmentMatches, nmatches, n, collect_fmatches)

        X, Hs, IDtoROW, last_matched_col = buildDesignMatrix(fragmentMatches, fragmentMisses, nmatches, nmisses, H_COLS, H_ROWS, H_VALS)

        weights = sparseNMF(Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)[:]

        scores = getDistanceMetrics(X, Hs, last_matched_col)
        
        unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches, nmatches,err_dist)

        Score!(scored_PSMs, unscored_PSMs, 
                length(spectra[:intensities][i]), 
                Float64(sum(spectra[:intensities][i])), 
                match_count/prec_count, scores[:scribe], scores[:city_block], scores[:matched_ratio], scores[:spectral_contrast], scores[:entropy_sim], weights, IDtoROW,
                scan_idx = i,
                min_spectral_contrast = min_spectral_contrast,
                min_frag_count = min_frag_count
                )
        
        #Reset 
        reset!(transitions, transition_idx)
        reset!(fragmentMatches, nmatches), reset!(fragmentMisses, nmisses)

    end

    if collect_frag_errs
        return DataFrame(scored_PSMs), all_matches
    else
        DataFrame(scored_PSMs)
    end

end

function reset!(fms::Vector{FragmentMatch{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        fms[i] = FragmentMatch{T}()
    end
end

function reset!(lf::Vector{LibraryFragment{T}}, last_non_empty::Int64) where {T<:AbstractFloat}
    for i in range(1, last_non_empty)
        lf[i] = LibraryFragment{T}()
    end
end


function collectFragErrs(all_fmatches::Vector{FragmentMatch{T}}, new_fmatches::Vector{FragmentMatch{T}}, nmatches::Int, n::Int, collect_fmatches::Bool) where {T<:AbstractFloat}
    if collect_fmatches
        for match in range(1, nmatches)
            if n < length(all_fmatches)
                all_fmatches[n] = new_fmatches[match]
                n += 1
            end
        end
    end
    return n
end

if collect_frag_errs
    for match in range(1, nmatches)
        if n < length(all_matches)
            all_matches[n] = fragmentMatches[match]
            n += 1
        end

    end
end
#=

function SearchRAW(
                    spectra::Arrow.Table, 
                    #ptable::PrecursorDatabase,
                    frag_index::FragmentIndex{Float32},
                    fragment_list::Vector{Vector{LibraryFragment{Float32}}},
                    ms_file_idx::UInt32,
                    interpolation::Any;
                    isolation_width::Float64 = 4.25,
                    precursor_tolerance::Float64 = 5.0,
                    fragment_tolerance::Float64 = 20.0,
                    frag_ppm_err::Float64 = 0.0,
                    topN::Int64 = 20,
                    min_frag_count::Int64 = 4,
                    min_matched_ratio::Float32 = Float32(0.8),
                    min_spectral_contrast::Float32 = Float32(0.65),
                    λ::Float32 = Float32(1e3),
                    γ::Float32 = zero(Float32),
                    scan_range::Tuple{Int64, Int64} = (0, 0), 
                    max_peaks::Int = 200, 
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0),
                    rt_tol::Float32 = Float32(30.0),
                    #fragment_match_ppm::U,
                    data_type::Type{T} = Float64,
                    sample_rate::Float64 = 1.0,
                    collect_frag_errs = false
                    ) where {T<:Real}
    
    scored_PSMs = makePSMsDict(XTandem(data_type))
    ms2, MS1, MS1_i = 0, 0, 0
    precs = Counter(UInt32, UInt8, Float32, length(fragment_list)) #Prec counter
    n = 1
    all_matches = Vector{FragmentMatch{Float32}}()
    if collect_frag_errs
        all_matches = [FragmentMatch{Float32}() for x in range(1, 30000)]
    end
    #for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    fragmentMatches = [FragmentMatch{Float32}() for x in range(1, 100000)]
    fragmentMisses = [FragmentMatch{Float32}() for x in range(1, 100000)]
    #for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))

        if (i%10000) == 0
            println(i)
        end
        if spectrum[:msOrder] == 1
            MS1 = spectrum[:masses]
            MS1_i = spectrum[:intensities]
            continue
        else
            ms2 += 1
        end

        if scan_range != (0, 0)
            i < first(scan_range) ? continue : nothing
            i > last(scan_range) ? continue : nothing
        end

        if first(rand(1)) > sample_rate
            continue
        end

        min_intensity = spectrum[:intensities][sortperm(spectrum[:intensities], rev = true)[min(max_peaks, length(spectrum[:intensities]))]]

        iRT = interpolation(min(spectrum[:retentionTime], 133.0))::Float64
        iRT_low = Float32(iRT - rt_tol)
        iRT_high = Float32(iRT + rt_tol)
        #reset!(precs)

        prec_count, match_count = searchScan!(precs,
                    frag_index, 
                    min_intensity, spectrum[:masses], spectrum[:intensities], spectrum[:precursorMZ], iRT_low, iRT_high,
                    Float32(fragment_tolerance), 
                    Float32(precursor_tolerance),
                    Float32(isolation_width),
                    min_frag_count = min_frag_count, 
                    min_ratio = Float32(min_matched_ratio),
                    topN = topN
                    )
        

        transitions = selectTransitions(fragment_list, precs, topN)

        reset!(precs)
        
        if length(transitions) == 0
            continue
        end
        
        nmatches, nmisses = matchPeaks(transitions, 
                                    fragmentMatches,
                                    fragmentMisses,
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    FragmentMatch{Float32},
                                    count_unmatched=true,
                                    #δs = zeros(T, (1,)),
                                    δs = [frag_ppm_err],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = min_intensity,
                                    ppm = fragment_tolerance
                                    )
        #=println("length(spectrum[:masses])", length(spectrum[:masses]))
        println("nmisses ", nmisses)
        println("nmatches ", nmatches)
        println("length(fragmentMatches) ", length(fragmentMatches))
        println("length(fragmentMisses) ", length(fragmentMisses))=#
        if nmatches < 2#iszero(length(fragmentMatches))
            for i in range(1, nmatches)
                fragmentMatches[i] = FragmentMatch{Float32}()
            end
    
            for i in range(1, nmisses)
                fragmentMisses[i] = FragmentMatch{Float32}()
            end
            continue
        end

        if collect_frag_errs
            for match in range(1, nmatches)
                if n < length(all_matches)
                    all_matches[n] = fragmentMatches[match]
                    n += 1
                end
    
            end
            #Eappend!(all_matches, fragmentMatches)
        end
 
        X, Hs, IDtoROW, last_matched_col = buildDesignMatrix(fragmentMatches, fragmentMisses, nmatches, nmisses)
        #println("size(Hs) ", size(Hs))
        #println("length(keys(IDtoROW))", length(keys(IDtoROW)))
        weights = sparseNMF(Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)[:]
        #return X, Hs, Hst, IDtoROW, weights

        #scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all, kt_pval = getDistanceMetrics(Hst, X, weights, matched_cols)
        scores = getDistanceMetrics(X, Hs, last_matched_col)
        
        #For progress and debugging. 

        unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches, nmatches)
        for i in range(1, nmatches)
            fragmentMatches[i] = FragmentMatch{Float32}()
        end

        for i in range(1, nmisses)
            fragmentMisses[i] = FragmentMatch{Float32}()
        end
        #times[:score] += @elapsed Score!(scored_PSMs, unscored_PSMs, 
        Score!(scored_PSMs, unscored_PSMs, 
                length(spectrum[:intensities]), 
                Float64(sum(spectrum[:intensities])), 
                match_count/prec_count, scores[:scribe], scores[:city_block], scores[:matched_ratio], scores[:spectral_contrast], scores[:entropy_sim], weights, IDtoROW,
                scan_idx = Int64(i),
                min_spectral_contrast = min_spectral_contrast,
                min_frag_count = min_frag_count
                )
    end

    if collect_frag_errs
        return DataFrame(scored_PSMs), all_matches
    else
        DataFrame(scored_PSMs)
    end
end

CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/psms_071923.csv",PSMs)

test_psms = SearchRAW(MS_TABLE, prosit_index_5ppm_15irt, frag_detailed, UInt32(1), linear_spline,
       min_frag_count = 4, 
       topN = 200, 
       fragment_tolerance = 15.6, 
       λ = Float32(1e4), 
       γ =Float32(1),
       max_peaks = 1000, 
       scan_range = (0, 300000), #101357 
       precursor_tolerance = 20.0,
       min_spectral_contrast =  Float32(0.65),
       rt_tol = Float32(15.0)
       )
2×18 DataFrame
 Row │ b_ladder  city_block  entropy     error       hyperscore  intensity_explained  matched_ratio  ms_file_idx  poisson   precursor_idx  scan_idx  scribe_score  spectral_contrast_all  spectral_contrast_matched  spectrum_peaks  total_ions  ⋯
     │ Int8      Float32     Float64     Float64     Float64     Float64              Float32        UInt32       Float64   UInt32         Int64     Float32       Float32                Float32                    UInt32          UInt32      ⋯
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │        6    -3.64194  -5.08692e8  0.00273705     94.9202            0.233232        5.26499             1  -83.0898        4211123    101357      15.9464                0.939872                   0.979832             642          32  ⋯
   2 │        1    -2.55285  -6.99821e7  0.00474243     38.5233            0.0327072       0.751578            1  -16.2589        6522764    101357       9.14359               0.817562                   0.918068             642          10
=#
   #=

    prosit_index_5ppm_5irt.fragment_bins[50893]
FragBin{Float32}(397.25577f0, 397.2571f0, 0x0000c6cd)
PSMs = DataFrame(CSV.File("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_071423.csv"))

 LibraryFragment{Float32}(397.25577f0, 0x01, true, 0x03, 0x01, 0.21178919f0, 0x03, 0x004041b3)
 LibraryFragment{Float32}(133.09012f0, 0x03, true, 0x03, 0x02, 0.0024077685f0, 0x03, 0x004041b3)
 LibraryFragment{Float32}(212.10297f0, 0x01, false, 0x03, 0x03, 0.19929862f0, 0x03, 0x004041b3)
 LibraryFragment{Float32}(468.2929f0, 0x01, true, 0x04, 0x04, 0.23923871f0, 0x03, 0x004041b3)

best_psms = combine(sdf -> sdf[argmin(sdf.RT_error),:], groupby(PSMs[(PSMs[:,:q_values].<=0.005).&(PSMs[:,:decoy].==false).&(PSMs[:,:RT_error].<1000.0),:], :precursor_idx))
best_psms = combine(sdf -> sdf[argmin(sdf.RT_error),:], groupby(PSMs[(PSMs[:,:matched_ratio].>10.0).&(PSMs[:,:decoy].==false).&(PSMs[:,:spectral_contrast_all].>0.97).&(PSMs[:,:q_values].<0.01),:], :precursor_idx))
plot(best_psms[:,:iRT], best_psms[:,:RT], seriestype=:scatter, xlabel = "iRT", ylabel = "RT")

ns1 = Splines2.ns_(collect(range(0.0, length=100, stop=130.0)),df=4,intercept=true);
X = ns1(best_psms[:,:RT]);
fit1 = lm(X,best_psms[:,:iRT]);
best_psms[:,:iRT_pred] =  GLM.predict(fit1, ns1(best_psms[:,:RT]))
plot(best_psms[:,:iRT], best_psms[:,:RT], seriestype=:scatter, xlabel = "iRT", ylabel = "RT")

#GLM.predict(fit1, ns1(collect(range(10.0, length = 100, stop = 130.0))))
plot!(GLM.predict(fit1, ns1(collect(range(0.0, length = 100, stop = 130.0)))),collect(range(10.0, length = 100, stop = 130.0)), xlabel = "iRT", ylabel = "RT")

sort_order = sortperm(best_psms[:,:iRT])
histogram(best_psms[sort_order[1000:1200],:RT])
#X, H, UNMATCHED, IDtoROW, fragmentMatches, fragmentMisses = 
test_psms = SearchRAW(MS_TABLE, prosit_index_5ppm_5irt, frag_detailed, UInt32(1), linear_spline,
       min_frag_count = 4, 
       topN = 200, 
       fragment_tolerance = 15.6, 
       λ = Float32(1e4), 
       γ =Float32(1),
       max_peaks = 1000, 
       scan_range = (101357, 101357), 
       precursor_tolerance = 20.0,
       min_spectral_contrast =  Float32(0.65)
       )


       using MultiKDE
using Distributions, Random, Plots
using MultiKDE
sort_order = sortperm(best_psms[:,:RT])
bins = [Int(x÷1) for x in range(1, length = 200, stop = length(sort_order))]
ys = []
xs = Float64[]
rts = Float64[]
for i in 1:(length(bins) - 1)
    observations = best_psms[sort_order[bins[i]:(bins[i + 1])],:iRT]
    x = Vector(LinRange(minimum(observations), maximum(observations), 50))
    kde = KDEUniv(ContinuousDim(), 3.0, observations, MultiKDE.gaussian)
    y = [MultiKDE.pdf(kde, _x, keep_all=false) for _x in x]
    push!(ys, y)
    push!(xs, x[argmax(y)])
    push!(rts, mean(best_psms[sort_order[bins[i]:(bins[i + 1])],:RT]))
    #push!(rts, best_psms[sort_order[bins[i]],:RT])
end
plot(best_psms[:,:RT], best_psms[:,:iRT], seriestype=:scatter, xlabel = "RT", ylabel = "iRT")
plot!(rts, xs)
linear_spline = LinearInterpolation(rts, xs, extrapolation_bc = Line() )
ns1 = Splines2.ns_(collect(range(0.0, length=30, stop=130.0)),df=10,intercept=false);
X = ns1(rts);
fit1 = lm(X,xs);
best_psms[:,:iRT_pred] =  linear_spline(best_psms[:,:RT])
plot!(best_psms[:,:RT], best_psms[:,:iRT_pred], xlabel = "RT", ylabel = "iRT")



highest = maximum([maximum(y) for y in ys])
plot(xs[1], ys[1], label=bws, fmt=:svg)
plot!(observations, [highest+0.05 for _ in 1:length(ys)], seriestype=:scatter, label="observations", size=(900, 450), legend=:outertopright)

# Simulation
# Simulation
bws = [1 3.0 5.0]
d = Normal(0, 1)
observations = best_psms[sort_order[1:200],:RT]
granularity_1d = 100
x = Vector(LinRange(0.0, 130.0, granularity_1d))
ys = []
for bw in bws
    kde = KDEUniv(ContinuousDim(), bw, observations, MultiKDE.gaussian)
    y = [MultiKDE.pdf(kde, _x, keep_all=false) for _x in x]
    push!(ys, y)
end

# Plot
highest = maximum([maximum(y) for y in ys])
plot(x, ys, label=bws, fmt=:svg)
plot!(observations, [highest+0.05 for _ in 1:length(ys)], seriestype=:scatter, label="observations", size=(900, 450), legend=:outertopright)
=#