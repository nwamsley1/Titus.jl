#using Dictionaries
#using Combinatorics
#include("src/precursor.jl")


"""
    digest(sequence::AbstractString; regex::Regex = r"[KR][^P|\$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1)

Given an amino acid sequence, `sequence`, and a regular expression that matches enzymatic cleavage sites, finds all peptide cleavage products from the
amino acid sequence. Can set a minimum and maximum length for cleaved peptides. Gets all peptides with 0:N missed cleavages. Default cleavage site
regex is for trypsin. Returns a list of strings that are the peptide cleavage products. 

- `sequence::AbstractString` -- An amino acid sequence
- `regex::Regex` -- Cleavage site regex. 
- `max_length::Int` -- Exclude peptides longer than this
- `min_length::Int` -- Exclude peptides shorter than this
- `missed_cleavages::Int` -- Get peptides with 0 to this many missed cleavages

"""
function digest(sequence::AbstractString; regex::Regex = r"[KR][^P|$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1)
    
    function addPeptide!(peptides::Vector{SubString{String}}, n::Int, sequence::AbstractString, site::Int, previous_sites::Vector{Int64}, min_length::Int, max_length::Int, missed_cleavages::Int)
        for i in 1:min(n, missed_cleavages + 1)
            previous_site = previous_sites[end - i + 1]
            if ((site - previous_site) >= min_length) && ((site - previous_site) <= max_length)
                    @inbounds push!(peptides, @view sequence[previous_site+1:site])
            end
        end

        @inbounds for i in 1:length(previous_sites)-1
            previous_sites[i] = previous_sites[i+1]
        end
        n += 1
        previous_sites[end] = site
        return n
    end

    #Iterator of cleavage sites in `sequence`
    cleavage_sites = eachmatch(regex, sequence, overlap = true)
    #in silico digested peptides from `sequence`
    peptides = Vector{SubString{String}}(undef, 0)
    #Positions in `sequence` of 1:(missed_cleavages + 1) cleavage sites
    previous_sites = zeros(Int64, missed_cleavages + 1)
    previous_sites[1] = 1
    #Number of cleavage sites encountered so far
    n = 1
    
    #Advance through each cleavage_site and see if there are any new peptides
    #that satisfy the max_length, min_length, and missed_cleavages thresholds. 
    for site in cleavage_sites
        n = addPeptide!(peptides, n, sequence, site.offset, previous_sites, min_length, max_length, missed_cleavages)
    end

    #Peptides containint the C-terminus are added here
    n = addPeptide!(peptides, n, sequence, length(sequence), previous_sites, min_length, max_length, missed_cleavages)
    return peptides 
end


#using FASTX
#using CodecZlib
#using Dictionaries
#file_path = "/Users/n.t.wamsley/RIS_temp/HAMAD_MAY23/mouse_SIL_List/UP000000589_10090.fasta.gz"
#function parseRule(identifier::String)
#    split(identifier, "|")[2]
#end

struct FastaEntry
    uniprot_id::String
    description::String
    sequence::String
    decoy::Bool
end

getID(fe::FastaEntry) = fe.uniprot_id
getDescription(fe::FastaEntry) = fe.description
getSeq(fe::FastaEntry) = fe.sequence
isDecoy(fe::FastaEntry) = fe.decoy


FastaEntry() = FastaEntry("","","", false)

function shufflefast(s::String)
 
    ss = sizeof(s)
    l = length(s) - 1

    v = Vector{Int}(undef, l)
    i = 1
    for j in 1:(l)
        v[j] = i
        i = nextind(s, i)
    end
    #println("TEST")
    p = pointer(s)
    u = Vector{UInt8}(undef, ss)
    k = 1
    for i in randperm(l)
        for j in v[i]:(i == l ? ss : v[i+1]-1)
            u[k] = unsafe_load(p, j)
            k += 1
        end
    end
    u[end] = unsafe_load(p, ss);
    return String(u)
end

"""
        parseFasta(fasta_path::String, parse_identifier::Function = x -> x)::Vector{FastaEntry}

    `    Given a fasta file, loads the identifier, description, and sequence into an in memory representation 
    `FastaEntry` for each entry. Returns Vector{FastaEntry}

    - `fasta_path::String` -- Path to a ".fasta" or ".fasta.gz"
    - `parse_identifier::Function ` -- Function takes a String argument and returns a String 
    `
"""
function parseFasta(fasta_path::String, parse_identifier::Function = x -> split(x,"|")[2])::Vector{FastaEntry}

    function getReader(fasta_path::String)
        if endswith(fasta_path, ".fasta.gz")
            return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
        elseif endswith(fasta_path, ".fasta")
            return FASTA.Reader(open(fasta_path))
        else
            throw(ErrorException("fasta_path \"$fasta_path\" did not end with `.fasta` or `.fasta.gz`"))
        end
    end

    #I/O for Fasta
    reader = getReader(fasta_path)

    #In memory representation of FastaFile
    fasta = Vector{FastaEntry}()
    @time begin
        for record in reader
                push!(fasta, 
                        FastaEntry(parse_identifier(FASTA.identifier(record)),
                                FASTA.description(record),
                                FASTA.sequence(record),
                                false)
                )
        end
    end

    return fasta
end

function digestFasta(fasta::Vector{FastaEntry}; regex::Regex = r"[KR][^P|$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1)
    peptides_fasta = Vector{FastaEntry}()
    #For eacbh protein in the FASTA
    for entry in fasta
        #Digest the protein into peptides
        for peptide in digest(getSeq(entry), regex = regex, max_length = max_length, min_length = min_length, missed_cleavages = missed_cleavages)
            #Make a target and decoy for each peptide
            for decoy in [false, true]
                if decoy
                    peptide = reverse(peptide[1:(end -1)])*peptide[end]#shufflefast(String(peptide))
                    #peptide = shufflefast(String(peptide))
                end
                push!(peptides_fasta, FastaEntry(getID(entry), 
                                                "",#getDescription(entry) This justs wastes time and memory here 
                                                peptide,
                                                decoy))
            end
        end
    end
    #Sort the peptides
    return sort!(peptides_fasta, by = x -> getSeq(x))
end

#digestFasta(parseFasta("/Users/n.t.wamsley/Projects/TEST_DATA/proteomes/UP000000589_10090.fasta.gz"))
#=
FASTA_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/"
fasta_paths = [joinpath(FASTA_DIR, file) for file in filter(file -> isfile(joinpath(FASTA_DIR, file)) && match(r".fasta.gz$", file) != nothing, readdir(FASTA_DIR))];
peptides_fasta = sort!(vcat([digestFasta(parseFasta(fasta)) for fasta in fasta_paths]...), by = x->getSeq(x))
precursor_table = PrecursorTable()
fixed_mods = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
var_mods = [(p=r"(M)", r="[MOx]")]
buildPrecursorTable!(precursor_table, peptides_fasta, fixed_mods, var_mods, 2)
const charge_facs = Float32[1, 0.9, 0.85, 0.8, 0.75]

function adjustNCE(NCE::T, default_charge::Integer, peptide_charge::Integer) where {T<:AbstractFloat}
    return NCE*(charge_facs[default_charge]/charge_facs[peptide_charge])
end

buildPrositCSV(precursor_table,  
                "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/humanYeastEcoli_NCE33_corrected_091623.csv",
                min_length = 7, 
                max_length = 30, 
                nce = 33.0, 
                dynamic_nce = true
                )

buildPrositCSV(precursor_table,  
                "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/humanYeastEcoli_NCE33_fixed_091623.csv",
                min_length = 7, 
                max_length = 30, 
                nce = 33.0, 
                dynamic_nce = false
                )



FASTA_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/"
fasta_paths = ["/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/UP000005640_9606_human.fasta.gz"]
peptides_fasta = sort!(vcat([digestFasta(parseFasta(fasta)) for fasta in fasta_paths]...), by = x->getSeq(x))
precursor_table = PrecursorTable()
fixed_mods = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
var_mods = [(p=r"(M)", r="[MOx]")]
buildPrecursorTable!(precursor_table, peptides_fasta, fixed_mods, var_mods, 2)
const charge_facs = Float32[1, 0.9, 0.85, 0.8, 0.75]

function adjustNCE(NCE::T, default_charge::Integer, peptide_charge::Integer) where {T<:AbstractFloat}
    return NCE*(charge_facs[default_charge]/charge_facs[peptide_charge])
end

buildPrositCSV(precursor_table,  
                "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/human_NCE33_corrected_091623.csv",
                min_length = 7, 
                max_length = 30, 
                nce = 33.0, 
                dynamic_nce = true
                )

buildPrositCSV(precursor_table,  
                "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/human_NCE33_fixed_091623.csv",
                min_length = 7, 
                max_length = 30, 
                nce = 33.0, 
                dynamic_nce = false
                )
=#