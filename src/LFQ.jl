using Statistics 

"""
    getS(peptides::AbstractVector{String}, peptides_dict::Dict{String, Int64}, experiments::AbstractVector{UInt32}, experiments_dict::Dict{UInt32, Int64}, abundance::AbstractVector{Union{T, Missing}}, M::Int, N::Int) where {T<:Real}

Get MxN matrix of intensities of each peptide in each experiment. Rows stand for experiments and columns for peptides. If the j'th peptide 
was not seen in the i'th experiment, then S[i, j] == missing. 

### Input 
- `peptides::AbstractVector{String}` -- MxN List of peptides for the protein. 
- `peptides_dict::Dict{String, Int64}` -- Maps N peptides by name to column numbers in S 
- `experiments::AbstractVector{UInt32}` -- MxN List of experiments 
- `experiments_dict::Dict{UInt32, Int64` -- Maps M experiment IDs by name to rows numbers in S 
- `abundance::AbstractVector{Union{T, Missing}}` -- MxN list of peptide abundances
- `M::Int` -- Number of experiments/rows in S
- `N::Int` -- Number of peptides/columns in S


### Examples 

"""
function getS(peptides::AbstractVector{String}, peptides_dict::Dict{String, Int64}, experiments::AbstractVector{UInt32}, experiments_dict::Dict{UInt32, Int64}, abundance::AbstractVector{Union{T, Missing}}, M::Int, N::Int) where {T<:Real}
    #Initialize
    S = Array{Union{Missing,T}}(undef, (M, N))
    for i in eachindex(peptides)
            if !ismissing(abundance[i]) #Abundance of the the peptide 
               if abundance[i] != 0.0
                    S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = abundance[i]
               else
                    S[peptides_dict[peptides[i]], experiments_dict[experiments[i]]] = missing
               end
            end
    end
    return S
end

"""
    getB(S::Matrix{Union{Missing, T}}, N::Int, M::Int) where {T<:Real}

Column vector of sum of median peptide ratios. See the following references. 

* Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
* https://rdrr.io/cran/iq/man/maxLFQ.html
* Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.

"""
function getB(S::Matrix{Union{Missing, T}}, N::Int, M::Int) where {T<:Real}
    B = zeros(N + 1)
    for i in 1:N
        for j in 1:N
            if j != i
                #Ratios of peptides commonly identified bewteen experiment i and j. 
                r_i_j = skipmissing(-log2.( @view(S[:,j]) ) .+ log2.(@view(S[:,i])))
                #r_i_j could be emptry if no peptides were commonly identified between experiment i and j. 
                if !isempty(r_i_j)
                    #median(r_i_j) is the median of all ratios between peptides commonly quantified between experiment
                    #i and experiment j. 
                    B[i] += median(r_i_j)
                end
            end
        end
    end
    B.*=2
    norm = 0.0

    #Relative abundances are correct without this step, but the absolute values should make sense. 
    for row in eachrow(S)
        norm += sum(log2.(skipmissing(row)))*(length(row)/(length(row) - sum(ismissing.(row))))
    end
    B[end] = norm/M #For scaling
    B
end

"""
    getA(N::Int)

Design matrix for maxLFQ. See the following

    * Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
    * https://rdrr.io/cran/iq/man/maxLFQ.html
    * Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.
    
"""
function getA(N::Int)
    A = ones(N+1, N+1)
    for i in 1:(N)
        for j in 1:N
            if i == j
                A[i, j] = 2*(N - 1)
            else
                A[i, j] = -2
            end
        end
    end
    A[end, end] = 0

    A
end

"""
    getProtAbundance(protein::String, peptides::AbstractVector{String}, experiments::AbstractVector{UInt32}, abundance::AbstractVector{Union{T, Missing}},
                            protein_out::Vector{String}, peptides_out::Vector{String}, experiments_out::Vector{UInt32}, log2_abundance_out::Vector{Float64}) where {T <: Real}

Estimates protein-level abundances from peptide level abundances using the MaxLFQ algorithm 
* Yu F, Haynes SE, Nesvizhskii AI. IonQuant Enables Accurate and Sensitive Label-Free Quantification With FDR-Controlled Match-Between-Runs. Mol Cell Proteomics. 2021;20:100077. doi: 10.1016/j.mcpro.2021.100077. Epub 2021 Apr 2. PMID: 33813065; PMCID: PMC8131922.
* https://rdrr.io/cran/iq/man/maxLFQ.html
* Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666.

### Input
`peptides`, `experiments`, and `abundance` each have length `N`
- `protein::String` -- Protien Name
- `peptides::AbstractVector{String}` -- Peptide names. Same length as `experiments` and `abundance`
- `experiments::AbstractVector{UInt32}` -- MS file id's. Same length as `peptides` and `abundance`
- `abundance::AbstractVector{Union{T, Missing}}` -- Abundances of the peptides. Same lenght as `peptides` and `experiments`
- `protein_out::Vector{String}` -- The protein name. Repeated `N` times for each experiment
- `peptides_out::Vector{String}` -- List of lists of detected peptides for each of the `N` experiments. [PEPA;PEPB PEPA;PEPC ...]
- `experiments_out::Vector{UInt32}` -- List of expeiments 
- `log2_abundance_out::Vector{Float64}` -- List of log2 estimated protein abundances for each experiment. 

### Examples. 
julia> using DataFrames

julia> prot = DataFrame(Dict(
           :peptides => ["A","A","A","B","B","B","C","C","C","D","D","D"],
           :protein => append!(split(repeat("A",9), ""), ["B","B","B"]),
           :file_idx => UInt32[1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],
           :abundance => [10, 20, 40, 1, 2, 4, 100, 200, missing, 1000, 2000, 3000],
       ))

12x4 DataFrame
 Row │ abundance  file_idx  peptides  protein   
     │ Int64?     UInt32    String    SubStrin… 
─────┼──────────────────────────────────────────
   1 │        10         1  A         A
   2 │        20         2  A         A
   3 │        40         3  A         A
   4 │         1         1  B         A
   5 │         2         2  B         A
   6 │         4         3  B         A
   7 │       100         1  C         A
   8 │       200         2  C         A
   9 │   missing         3  C         A
  10 │      1000         1  D         B
  11 │      2000         2  D         B
  12 │      3000         3  D         B

julia> function testLFQ(prot)
           out = Dict(
               :protein => String[],
               :peptides => String[],
               :log2_abundance => Float64[],
               :experiments => UInt32[],
           )

           for (protein, data) in pairs(groupby(prot, :protein))
               getProtAbundance(string(protein[:protein]), 
                                   collect(data[!,:peptides]), 
                                   collect(data[!,:file_idx]), 
                                   collect(data[!,:abundance]),
                                   out[:protein],
                                   out[:peptides],
                                   out[:experiments],
                                   out[:log2_abundance]
                               )
           end
           out
       end
testLFQ (generic function with 1 method)

julia> out = testLFQ(prot)
Dict{Symbol, Vector} with 4 entries:
  :protein        => ["A", "A", "A", "B", "B", "B"]
  :peptides       => ["A;B;C", "A;B;C", "A;B", "D", "D", "D"]
  :log2_abundance => [3.15526, 4.15526, 5.15526, 9.96578, 10.9658, 11.5507]
  :experiments    => UInt32[0x00000001, 0x00000002, 0x00000003, 0x00000001, 0x00000002, 0x00000003]

"""
function getProtAbundance(protein::String, peptides::AbstractVector{String}, experiments::AbstractVector{UInt32}, abundance::AbstractVector{Union{T, Missing}},
                          protein_out::Vector{String}, peptides_out::Vector{String}, experiments_out::Vector{UInt32}, log2_abundance_out::Vector{Float64}) where {T <: Real}

    unique_experiments = unique(experiments)
    unique_peptides = unique(peptides)

    N = length(unique_experiments)
    M = length(unique_peptides)

    peptides_dict = Dict(zip(unique_peptides, 1:M))
    experiments_dict = Dict(zip(unique_experiments, 1:N))

    #Appends the results to the inputs `protein_out`, `peptides_out`, `experiments_out` and `log2_abundance_out`
    function appendResults!(protein::String, peptides::Vector{String}, log2_abundances::Vector{T}, protein_out::Vector{String}, peptides_out::Vector{String}, log2_abundance_out::Vector{T}, S::Matrix{Union{Missing, U}}) where {T,U<:Real}
        
        function appendPeptides!(peptides_out::Vector{String}, peptides::Vector{String}, S::Matrix{Union{Missing,U}})
            #Each column of S corresponds to and experiment and each row corresponds to a peptide
            #Need to get each non-missing peptide for each experiment. Concatenate the non-missing peptides
            #For and experiment with a semi-colon. See example in the getProtAbundance docstring 
            for j in eachindex(eachcol(S))
                peps = String[]
                for i in eachindex(@view(S[:,j]))
                    if !ismissing(S[i,j])
                        push!(peps, peptides[i])
                    end
                end
                push!(peptides_out, join(peps, ";"))
            end
        end
        append!(log2_abundance_out, log2_abundances)
        append!(experiments_out, unique_experiments)
        append!(protein_out, [protein for x in 1:N])
        appendPeptides!(peptides_out, peptides, S)
    end

    #ixj matrix where rows are for experiments and columns are for peptides. Each entry is the abundance of the peptide
    #in the given experiment, or missing if peptide j was not seen in experiment i. 
    S = allowmissing(getS(peptides, peptides_dict, experiments, experiments_dict, abundance, M, N))

    #Column vector. The response matrix in Ax=B
    B = getB(S, N, M)

    #Design matrix. See references in docstring. 
    A = getA(N)

    #Solve linear system to get log-2 abundances 
    log2_abundances = (A\B)[1:(end - 1)]
    appendResults!(protein, unique_peptides, log2_abundances, protein_out, peptides_out, log2_abundance_out, S)

end
