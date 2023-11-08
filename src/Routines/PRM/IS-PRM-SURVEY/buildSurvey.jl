#julia ./src/Routines/PRM/IS-PRM-SURVEY/buildSurvey.jl ./data/example_config/BUILD-SURVEY-TEST.json /Users/n.t.wamsley/Desktop/PROTEINPEPTIDE.txt /Users/n.t.wamsley/Desktop

using JSON, PrettyPrinting, PDFmerger, ArgParse, Arrow, DataFrames, Tables, CSV
#=
example params_json file

{
        "fixed_mods":[
                     ["C","C[Carb]"],
                     ["K$","K[Hlys]"],
                     ["R$","R[Harg]"]
        ],
        "modification_masses":
        {
        "Carb":57.021464,
        "Harg":10.008269,
        "Hlys":8.014199
        },
}
=#
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "params_json"
            help = "Path to a .json file with the parameters"
            required = true
        "precursor_list"
            help = "Path to a tab delimited table of precursors"
            required = true
        "out_dir"
            help = "Path to a folder where results will be saved"
            required = true
    end

    return parse_args(s)
end

ARGS = parse_commandline()

params = JSON.parse(read(ARGS["params_json"], String))
OUT_DIR = ARGS["out_dir"]
PRECURSOR_LIST_PATH = ARGS["precursor_list"]

out_folder = joinpath(OUT_DIR, "SurveyMethod")
if !isdir(out_folder)
    mkpath(out_folder)
end

function parse_mods(fixed_mods)
    fixed_mods_parsed = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
    for mod in fixed_mods
        push!(fixed_mods_parsed, (p=Regex(mod[1]), r = mod[2]))
    end
    return fixed_mods_parsed
end

fixed_mods = parse_mods(params["fixed_mods"])
modification_masses = Dict{String, Float64}(k => Float64(v) for (k, v) in params["modification_masses"])
charges = params["charges"]

using Arrow, DataFrames, Tables
using Plots

struct SurveyPrecursor{T<:AbstractFloat,U<:Unsigned}
    Compound::String
    mz::T
    z::U
    intensity_threshold::T
end

[include(joinpath(pwd(), "src", jl_file)) for jl_file in ["precursor.jl","binaryRangeQuery.jl","matchpeaks.jl","getPrecursors.jl","SearchRAW.jl","applyMods.jl"]]
#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","PRM", jl_file)) for jl_file in ["precursorChromatogram.jl","getMS1PeakHeights.jl"]]
#Files needed for PSM scoring
[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","FastXTandem.jl"]]
#Files needed for IS-PRM-SURVEY routine
[include(joinpath(pwd(), "src", "Routines","PRM","IS-PRM-SURVEY", jl_file)) for jl_file in ["initTransitions.jl","selectTransitions.jl","getBestTransitions.jl","buildPrecursorTable.jl","getBestPSMs.jl","writeTables.jl","plotPRM.jl"]]


#=
fixed_mods = [(p=r"C", r="C[Carb]"), (p=r"(K\$)", r="K[Hlys]"), (p=r"(R\$)", r="R[Harg]")]
mods_dict = Dict{String, Float64}("Carb" =>57.021464,
            "Harg" =>0.008269,
            "Hlys" =>8.014199
            )
survey_df = parsePrecursorList("/Users/n.t.wamsley/Desktop/PROTEINPEPTIDE.txt", fixed_mods, mods_dict, UInt8[2, 3])
=#
function parsePrecursorList(f_path::String,
                            fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                            mods_dict::Dict{String, <:AbstractFloat},
                            charges::Vector{<:Unsigned})
    survey_precursors = Vector{SurveyPrecursor{Float64, UInt8}}()
    open(f_path) do f
        for (row, protein_peptide) in enumerate(eachline(f))
            line = map(string, split(protein_peptide, "\t"))
            protein, peptide = line[1:2];
            peptide = fixedMods(peptide, fixed_mods);
            for charge in charges
                prec = Precursor(peptide, 
                                    mods_dict = mods_dict,
                                    charge = charge,
                                    isotope = UInt8(0))
                push!(survey_precursors,
                SurveyPrecursor(protein*"_"*peptide,
                                getMZ(prec),
                                charge,
                                1e4)
                )
            end
        end
    end
    survey_df = DataFrame(survey_precursors)
    survey_df[!,:Formula] .= ""
    survey_df[!,:Adduct] .= ""
    rename!(survey_df, :mz => Symbol("m/z"))
    return survey_df[!,[:Compound,:Formula,:Adduct,Symbol("m/z"),:z,:intensity_threshold]]
end

survey_df = parsePrecursorList(PRECURSOR_LIST_PATH, 
                                fixed_mods, 
                                modification_masses, 
                                UInt8.(charges)
                                );
CSV.write(joinpath(out_folder,"SurveyMethod.csv"), survey_df);