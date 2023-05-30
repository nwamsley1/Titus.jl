var documenterSearchIndex = {"docs":
[{"location":"introduction/","page":"Titus","title":"Titus","text":"CurrentModule = Titus","category":"page"},{"location":"introduction/#Titus","page":"Titus","title":"Titus","text":"","category":"section"},{"location":"introduction/","page":"Titus","title":"Titus","text":"Documentation for Titus.","category":"page"},{"location":"introduction/#A-Header","page":"Titus","title":"A Header","text":"","category":"section"},{"location":"introduction/#Another-Header","page":"Titus","title":"Another Header","text":"","category":"section"},{"location":"README/#Aims","page":"Introduction","title":"Aims","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"Titus aims to simplify preparation and analysis of IS-PRM experiments [1]. Supports analysis of survey runs to characterize internal standard peptides, preparation of methods, peak area ratio estimation and protein-level quantitation. Generates high-quality, informative chromatogram and spectra plots in multi-page pdf's, which are ctrl+f searchable in any standard pdf viewer. Titus is fast, supports multi-threading, and enables analysis of >100 experiments possible in ~1 min including compilation time. Chromatogram plot generation make take a bit longer. Still a work-in-progress and feedback welcome. ","category":"page"},{"location":"README/#Features","page":"Introduction","title":"Features","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"Cross platform (tested on MacOS and ubuntu)\nAccepts raw ms data in the Apache Arrow format. See the following for cross-platform conversion of Thermo .raw files to the arrow format https://github.com/nwamsley1/ThermoRawFileToParquetConverter.\nAnalysis of survey methods. Given a list table of protein-peptide pairs, identifies the best charge state for each precursor (by XTandem hyperscore), the best transitions, and the MS1 peak height. If the survey analyses are split accross multiple experiments, these can be analyzed at once and combined. In addition, can run survey analyses at multiple collision energies/FAIMS CV's to identify the optimum for each analyte. Output is given in a format freindly to XCalibur method editor for Thermo Tribrid instruments.\nSupports variable and fixed modifications defined by regular expressions and includes examples. \nEstimates peak area ratios using an MM-Estimator (https://github.com/getzze/RobustModels.jl/blob/main/docs/make.jl). Enables accurate par estimation in the pressence of noisy or interfered transitions. High uncertainty in estimation can be used as grounds for exclusion. \nSummarizaiton of peptide-level quantitation to protein-level quantitation using the MaxLFQ Algorithm without normalization [2].\nGenerates a multi-page pdf for each experiment file including chromatogram plots. ","category":"page"},{"location":"README/","page":"Introduction","title":"Introduction","text":"(Image: alt text)","category":"page"},{"location":"README/#Future-Additions/In-Progress","page":"Introduction","title":"Future Additions/In-Progress","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"Work in progress \nMore complete documentation. At present only basic usage examples are given, and docstrings for most methods are listed in the CI docs but in no particular order. https://documentation.divio.com/\nRobust benchmarking of PAR estimation using robust linear estimators agasint against previously used methods, such as top-N [3] or exclusion of transitions based on cosine-similarity metrics [4].\nIf standard curves are available for absolute quantitation, the ability to convert PAR estimates to abundances\nBetter plot annotations\nCurrently not available as an Julia package. Usability is akward and could use improvement. Not all parameters can be user-defined yet. \nSome performance issues in \"precursors.jl\" that should be reasonably easy to resolve. \nFuture support of standard PRM and DIA using prosit libraries and spectral devonvolution?\nLots of others. Suggest your own....","category":"page"},{"location":"README/#Usage","page":"Introduction","title":"Usage","text":"","category":"section"},{"location":"README/#Survey-Run-Analyses","page":"Introduction","title":"Survey Run Analyses","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"julia ./src/Routines/PRM/IS-PRM_SURVEY/routine.jl ./data/test.json ./data/parquet/ ./data/NRF2_SIL.txt","category":"page"},{"location":"README/","page":"Introduction","title":"Introduction","text":"|Name                |Default| Short        |Description                    |  |––––––––––|–––-|––––––-|––––––––––|  |paramsjson||mandatory|Path to a .json file with the parameters (see Configuration)  |datadir||mandatory|\"Path to a folder with .arrow MS data tables\"  |precursorlist||mandatory|\"Path to a tab delimited table of precursors\"  |–makeplots|true|-p|\"Whether to make plots. Defaults to true\"  |–print_params|false|-s|\"Whether to print the parameters from the json. Defaults to false\"","category":"page"},{"location":"README/#Precursor-List-Example","page":"Introduction","title":"Precursor List Example","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"ABCB6\tYYNAESYEVER\nABCB6\tIDGQDISQVTQASLR\nABCB6\tALNVLVPIFYR\nABHD4\tYVSLPNQNK\nADD2\tVNVADEVQR\n.\n.\n.","category":"page"},{"location":"README/#IS-PRM-Analysis","page":"Introduction","title":"IS-PRM Analysis","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"julia --threads 24 ./src/Routines/PRM/IS-PRM/routine.jl ./data/IS-PRM_TEST.json ./data/parquet ./data/parquet/transition_list.csv","category":"page"},{"location":"README/","page":"Introduction","title":"Introduction","text":"|Name                |Default| Short        |Description                    |  |––––––––––|–––-|––––––-|––––––––––|  |paramsjson||mandatory|Path to a .json file with the parameters (see Configuration)  |datadir||mandatory|\"Path to a folder with .arrow MS data tables\"  |transitionlist||mandatory|\"Path to a tab delimited table of transitions\"  |–makeplots|true|-p|\"Whether to make plots. Defaults to true\"  |–print_params|false|-s|\"Whether to print the parameters from the json. Defaults to false\"","category":"page"},{"location":"README/#Transition-List-Example","page":"Introduction","title":"Transition List Example","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"protein_name,sequence,precursor_charge,precursor_isotope,transition_names\nABCB6,YYNAESYEVER[Harg],2,0,y10+1;y7+1;y8+1;y9+1;y6+1\nABCB6,ALNVLVPIFYR[Harg],2,0,y6+1;y8+1;y7+1;b4+1;y9+1\nABHD4,YVSLPNQNK[Hlys],2,0,b4+1;y7+1;y5+2;y3+1;y5+1\n.\n.\n.","category":"page"},{"location":"README/#Example-Outputs","page":"Introduction","title":"Example Outputs","text":"","category":"section"},{"location":"README/#Configuration-files","page":"Introduction","title":"Configuration files","text":"","category":"section"},{"location":"README/#IS-PRM-Survey","page":"Introduction","title":"IS-PRM-Survey","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"{\n    \"right_precursor_tolerance\": 0.001,\n    \"left_precursor_tolerance\": 0.001,\n    \"precursor_rt_tolerance\": 0.3,\n    \"b_ladder_start\": 3,\n    \"y_ladder_start\": 3,\n    \"precursor_charges\": [2, 3, 4],\n    \"precursor_isotopes\": [0],\n    \"transition_charges\": [1, 2],\n    \"transition_isotopes\": [0],\n    \"fragment_match_ppm\": 40,\n    \"minimum_fragment_count\": 5,\n    \"fragments_to_select\": 5,\n    \"precursor_rt_window\": 0.3,\n    \"max_variable_mods\": 2,\n    \"fixed_mods\":[\n                     [\"C\",\"C[Carb]\"],\n                     [\"K$\",\"K[Hlys]\"],\n                     [\"R$\",\"R[Harg]\"]\n    ],\n    \"variable_mods\":\n    [],\n    \"modification_masses\":\n        {\n        \"Carb\":57.021464,\n        \"Harg\":10.008269,\n        \"Hlys\":8.014199\n        },\n    \"ms_file_conditions\":\n        {\n            \"_35NCE_\":\"35NCE\",\n            \"_40NCE_\":\"40NCE\",\n            \"GAPDH\":\"GAPDH\"\n        }\n}","category":"page"},{"location":"README/#IS-PRM","page":"Introduction","title":"IS-PRM","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"{\n    \"right_precursor_tolerance\": 0.001,\n    \"left_precursor_tolerance\": 0.001,\n    \"precursor_rt_tolerance\": 0.3,\n    \"b_ladder_start\": 3,\n    \"y_ladder_start\": 3,\n    \"precursor_charges\": [2, 3, 4],\n    \"precursor_isotopes\": [0],\n    \"transition_charges\": [1, 2],\n    \"transition_isotopes\": [0],\n    \"fragment_match_ppm\": 40,\n    \"minimum_fragment_count\": 5,\n    \"fragments_to_select\": 5,\n    \"precursor_rt_window\": 0.3,\n    \"max_variable_mods\": 2,\n    \"fixed_mods\":[\n                     [\"C\",\"C[Carb]\"],\n                     [\"K$\",\"K[Hlys]\"],\n                     [\"R$\",\"R[Harg]\"]\n    ],\n    \"variable_mods\":\n    [],\n    \"modification_masses\":\n        {\n        \"Carb\":57.021464,\n        \"Harg\":10.008269,\n        \"Hlys\":8.014199\n        },\n    \"ms_file_conditions\":\n        {\n            \"_35NCE_\":\"35NCE\",\n            \"_40NCE_\":\"40NCE\",\n            \"GAPDH\":\"GAPDH\"\n        }\n}","category":"page"},{"location":"README/#References","page":"Introduction","title":"References","text":"","category":"section"},{"location":"README/","page":"Introduction","title":"Introduction","text":"<a id=\"1\">[1]</a>  Gallien S, Kim SY, Domon B. Large-Scale Targeted Proteomics Using Internal Standard Triggered-Parallel Reaction Monitoring (IS-PRM). Mol Cell Proteomics. 2015 Jun;14(6):1630-44. doi: 10.1074/mcp.O114.043968. Epub 2015 Mar 9. PMID: 25755295; PMCID: PMC4458725. <br> <a id=\"1\">[2]</a>  Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666 <br> <a id=\"1\">[3]</a>  Stopfer LE, Flower CT, Gajadhar AS, Patel B, Gallien S, Lopez-Ferrer D, White FM. High-Density, Targeted Monitoring of Tyrosine Phosphorylation Reveals Activated Signaling Networks in Human Tumors. Cancer Res. 2021 May 1;81(9):2495-2509. doi: 10.1158/0008-5472.CAN-20-3804. Epub 2021 Jan 28. PMID: 33509940; PMCID: PMC8137532. <br> <a id=\"1\">[4]</a>  Wamsley et al. Targeted proteomic quantitation of NRF2 signaling and predictive biomarkers in HNSCC. https://doi.org/10.1101/2023.03.13.532474 ","category":"page"},{"location":"","page":"Index","title":"Index","text":"CurrentModule = Titus","category":"page"},{"location":"#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"","page":"Index","title":"Index","text":"","category":"page"},{"location":"","page":"Index","title":"Index","text":"Modules = [Titus]","category":"page"},{"location":"#Titus.FragmentMatch","page":"Index","title":"Titus.FragmentMatch","text":"FragmentMatch\n\nType that represents a match between a fragment ion and a mass spectrum peak\n\nFields\n\ntransition::Transition – Represents a fragment ion\nintensity::Float32 – Intensity of the matching peak\nmatch_mz::Float32 – M/Z of the matching empirical peak.    NOT the transition mass and may differ from getMZ(transition) by some error\ncount::UInt8 – Number of matches (may want to count the number of matches if the spectrum is researched at    different mass offsets, such as in a cross correlation score)\npeak_int::Int64 – Index of the matching peak in the mass spectrum. \n\nExamples\n\n`FragmentMatch(transition::Transition, intensity::Float32, mass::Float32, count::UInt8, peak_ind::Int64) –    default internal constructor\nFragmentMatch() – constructor for null/empty precursor\n\nGetterMethods\n\ngetMZ(f::FragmentMatch) = getMZ(f.transition)\ngetLow(f::FragmentMatch) = getLow(f.transition)\ngetHigh(f::FragmentMatch) = getHigh(f.transition)\ngetPrecID(f::FragmentMatch) = getPrecID(f.transition)\ngetCharge(f::FragmentMatch) = getCharge(f.transition)\ngetIsotope(f::FragmentMatch) = getIsotope(f.transition)\ngetIonType(f::FragmentMatch) = getIonType(f.transition)\ngetInd(f::FragmentMatch) = getInd(f.transition)\n\n\n\n\n\n","category":"type"},{"location":"#Titus.Ion","page":"Index","title":"Titus.Ion","text":"Ion\n\nAbstract type that represents an ion \n\nTypes that inherit from Ion should implement the following. \n\ngetMZFeature(ion::Ion) = ion.mz\ngetMZ(ion::Ion) = getMZ(getMZFeature(ion))\ngetLow(ion::Ion) = getLow(getMZFeature(ion))\ngetHigh(ion::Ion) = getHigh(getMZFeature(ion))\ngetPrecID(ion::Ion) = ion.prec_id\ngetCharge(ion::Ion) = ion.charge\ngetIsotope(ion::Ion) = ion.isotope\n\n\n\n\n\n","category":"type"},{"location":"#Titus.MzFeature","page":"Index","title":"Titus.MzFeature","text":"MzFeature\n\nType that represents an m/z ratio with a ppm tolerance\n\nFields\n\nmono::Float32 – MZ with upper and lower bounds given a ppm tolerance\nlow::Float32 – Identifier of the precursor ion (parent ion of the fragment/transition)\nhigh::Float32 – Type of transition. For example 'b' for b ion or 'y' for y ion\n\nExamples\n\nMzFeature(mono::Float32; ppm::Float32 = Float32(20)) – default internal constructor\nMzFeature() = MzFeature(Float32(0.0)) – constructor for a default/placeholder MzFeature\n\nGetterMethods\n\ngetMZ(mzfeature::MzFeature) = mzfeature.mono\ngetLow(mzfeature::MzFeature) = mzfeature.low\ngetHigh(mzfeature::MzFeature) = mzfeature.high\n\n\n\n\n\n","category":"type"},{"location":"#Titus.Precursor","page":"Index","title":"Titus.Precursor","text":"Precursor <: Ion\n\nType that represents a precursor (A peptide parent ion)\n\nFields\n\nresidues::Vector{Residue} – List of amino acid residues of the precursor in their appropriate order\nmz::MzFeature – MZ with upper and lower bounds given a ppm tolerance\nprec_id::UInt32 – Identifier of the precursor ion (parent ion of the fragment/transition)\ncharge::UInt8 – Charge of the fragment ion\nisotope::UInt8 – Difference in number of neutrons from the monoisotopic fragment \npep_id::UInt32 – Identifier of the peptide from which precursor is derived\n\nExamples\n\nPrecursor(residues::Vector{Residue}, mz::Float32, charge::UInt8,             isotope::UInt8,             pep_id::UInt32,            prec_id::UInt32;             ppm = Float32(20)) – default internal constructor\nPrecursor() – constructor for null/empty precursor\nPrecursor(residues::Vector{Residue}, charge::UInt8,             isotope::UInt8 = UInt8(0),             pep_id::UInt32 = UInt32(0),            prec_id::UInt32 = UInt32(0) – Constructor that calculates mz without having to supply it\nPrecursor(sequence::String; mods_dict::Dict{String, Float32} = Dict{String, Float32}(), charge::UInt8 = UInt8(2),             isotope::UInt8 = UInt8(0),             pep_id::UInt32 = UInt32(0),            prec_id::UInt32 = UInt32(0)) – Constructor that accepts a string representation of a peptide\n\nGetterMethods\n\ngetResidues(precursor::Precursor) = precursor.residues\n\n\n\n\n\n","category":"type"},{"location":"#Titus.Precursor-Tuple{}","page":"Index","title":"Titus.Precursor","text":"Precursor()\nConstructor for an \"empty\" or \"default\" precursor\n\n\n\n\n\n","category":"method"},{"location":"#Titus.Precursor-Union{Tuple{String}, Tuple{T}} where T<:AbstractFloat","page":"Index","title":"Titus.Precursor","text":"Precursor(sequence::String, mods_dict::Dict{String, Float32}, charge::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0))\n\nAlternate constructor for the Precursor struct. Can accept a string representation of a peptide and a mods_dict  and convert to residues Array{Residue, 1}.  (link to Precusor)\n\n\n\n\n\n","category":"method"},{"location":"#Titus.Precursor-Union{Tuple{T}, Tuple{Array{Residue{T}, 1}, UInt8}, Tuple{Array{Residue{T}, 1}, UInt8, UInt8}, Tuple{Array{Residue{T}, 1}, UInt8, UInt8, UInt32}, Tuple{Array{Residue{T}, 1}, UInt8, UInt8, UInt32, UInt32}} where T<:AbstractFloat","page":"Index","title":"Titus.Precursor","text":"Precursor(residues::Array{Residue, 1}, charge::Int32, isotope::Int32 = Int32(0), prec_id::Int32 = Int32(0), pep_id::Int32 = Int32(0))\n\nConstructor for the Precursor struct. Given a list of amino acid residues, a charge, and an isotope state, makes a precursor object with the correct mz.  (link to Precursor)\n\n\n\n\n\n","category":"method"},{"location":"#Titus.Residue","page":"Index","title":"Titus.Residue","text":"Residue\n\nType that represents a (potentially modified) amino acid within a peptide\n\nFields\n\nmass:Float32 – mass of the amino acid\n\nExamples\n\nResidue(aa::AA) – default constructor\nResidue(aa::Char) = Residue(AA(aa))\nResidue(aa::AA, mod::Mod) = Residue(getMass(mod)+getMass(aa))\nResidue(residue::String, mods_dict::Dict{String, Float32}) = Residue(AA(residue[1]), Mod(residue, mods_dict))\nResidue(residue::String) = Residue(residue, Dict{String, Float32}())\nResidue(residue::String, mod_mass::Float32) = Residue(getMass(AA(residue[1])) + mod_mass)\nResidue(residue::Char, mod_mass::Float32) = Residue(getMass(AA(residue)) + mod_mass)\n\nGetter methods\n\ngetMass(residue::Residue) = residue.mass\n\nSee Also\n\ngetResidues(sequence::String,               mods_dict::Dict{String, Float32} = default_mods) – Gets a vector of residues given a string  representation of an amino acid sequence\n\n\n\n\n\n","category":"type"},{"location":"#Titus.Transition","page":"Index","title":"Titus.Transition","text":"Transition <: Ion\n\nType that represents transition (fragment ion of a peptide)\n\nFields\n\nmz::MzFeature – MZ with upper and lower bounds given a ppm tolerance\nprec_id::UInt32 – Identifier of the precursor ion (parent ion of the fragment/transition)\nion_type::Char  – Type of transition. For example 'b' for b ion or 'y' for y ion\nind::UInt8 – Position of fragment ion with reference to parent ion. (A b5+2 ion should have an ind equal to 5)\ncharge::UInt8 – Charge of the fragment ion\nisotope::UInt8 – Difference in number of neutrons from the monoisotopic fragment \n\nExamples\n\nTransition(frag_mz::Float32, prec_id::UInt32, ion_type::Char, ind::UInt8,              charge::UInt8,              isotope::UInt8;              ppm = Float32(20)) – default internal constructor\nTransition(residues::Vector{Residue}; ion_type::Char = 'y', charge::UInt8 = UInt8(1),              ind::UInt8 = UInt8(length(residues)),              isotope::UInt8 = UInt8(0),              prec_id::UInt32 = UInt32(0)) – constructor that calculates the appropriate mz\n\nGetterMethods\n\ngetIonType(transition::Transition) = transition.ion_type\ngetInd(transition::Transition) = transition.ind\n\n\n\n\n\n","category":"type"},{"location":"#Titus.Transition-Union{Tuple{Array{Residue{T}, 1}}, Tuple{T}} where T<:AbstractFloat","page":"Index","title":"Titus.Transition","text":"Transition(residues::Vector{Residue}; ion_type::Char = 'y', charge::UInt8 = UInt8(1), ind::UInt8 = UInt8(length(residues)), isotope::UInt8 = UInt8(0), prec_id::UInt32 = UInt32(0))\n\nConstructor for the `Transition` struct. Given a list of amino acid residues, ion type (b or y), charge state, and isotopic state makes a transition with the correct mz.\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getBIonModifier-Union{Tuple{UInt8}, Tuple{T}, Tuple{UInt8, T}} where T<:AbstractFloat","page":"Index","title":"Titus.getBIonModifier","text":"getBIonModifier(charge::UInt8)\n\nMass modification that needs to be added to b-ions\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getFragIons-Union{Tuple{Array{Residue{T}, 1}}, Tuple{T}} where T<:AbstractFloat","page":"Index","title":"Titus.getFragIons","text":"getFragIons(residues::Vector{Residue}; charge::UInt8 = UInt8(1), isotope::UInt8 = UInt8(0), y_start::Int = 3, b_start::Int = 3)\n\nUses calls to getIonSeries to concatenate both the b and y ion series together in a single Vector{Float32}\n\nInput\n\nresidues::Vector{Residue}: – List of amino acid residues in the peptide ion\ncharge::UInt8 – Charge of the fragment ions\nb_start::Int=3  – Index of first ion the the b-ion series to compute. \ny_start::Int=3  – Index of first ion the the y-ion series to compute. \nisotope::UInt8=UInt8(0) – Diference in the number of isotopes from the monoisotopic ion. \n\nOutput\n\nA Vector{Float32} wich each m/z in the ion series\n\nNotes\n\nThe modifier argument ought to depend on the kind of ion. For a 'b' ion series PROTON*charge is appropriate,\n\nbut for a 'y' ion series, PROTON*charge + H2O would be appropriate. \n\nWill not allow the index of the ion to be equal to or less than the charge. For example,\n\nb2+2 ions could only be calculated in error and are therefore excluded even if start  is set to 2. \n\nIf start exceeds length(residues)-1, then only the N-1 ion is calculated, that is,\n\nthe highest mass ion in the series. \n\nExamples\n\n#Gets the b3+1-b6+1 and y3+1-y6+1 ions\n\njulia> getFragIons(reverse(getResidues(\"PEPTIDE\")), b_start = 3, y_start = 3)\n8-element Vector{Float32}:\n 358.16083\n 459.2085\n 556.2612\n 685.30383\n 342.16595\n 443.21365\n 556.29767\n 671.3246\n\nSee Also\n\nAlternate convience methods\n\ngetFragIons(residues::Vector{Residue}; charge::UInt8 = UInt8(1),             isotope::UInt8 = UInt8(0),              y_start::Int = 3,              b_start::Int = 3) - Default method\ngetFragIons(precursor::Precursor; charge::UInt8 = UInt8(1),              isotope::UInt8 = UInt8(0),              y_start::Int = 3,              b_start::Int = 3) - Can supply a Precursor rather than a  Vector{Residues} input\ngetFragIons(precursor::Precursor,charges::Vector{UInt8},              isotopes::Vector{UInt8};              y_start::Int = 3,              b_start::Int = 3) - Gets b and y ion seriers for multiple charge and isotopic states\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getIonMZ-Union{Tuple{T}, Tuple{Array{Residue{T}, 1}, Char, UInt8}} where T<:AbstractFloat","page":"Index","title":"Titus.getIonMZ","text":"getIonMZ(residues::Vector{Residue}, ion_type::Char, charge::UInt8; isotope::UInt8 = UInt8(0))\n\nAlternate getIonMZ method that chooses the correct mass modifier for 'b', 'y', and 'p' ions respectively. \n\nInput\n\n- `residues::Vector{Residue}`: -- List of amino acid residues in the peptide ion\n- `ion_type::Char` -- Type of fragment ion. Currently supports, 'b', 'y', and 'p'. \n- `charge::UInt8` -- Charge of the ion\n- `isotope::UInt8=UInt8(0)` -- Diference in the number of isotopes from the monoisotopic ion.\n\nNotes\n\nSee main method\n    getIonMZ(residues::Vector{Residue}, charge::UInt8; modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))\nfor more details\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getIonMZ-Union{Tuple{T}, Tuple{Array{Residue{T}, 1}, UInt8}} where T<:AbstractFloat","page":"Index","title":"Titus.getIonMZ","text":"getIonMZ(residues::Vector{Residue}, charge::UInt8; modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))::Float32\n\nGet the mz ratio of an ion\n\nInput\n\nresidues::Vector{Residue}: – List of amino acid residues in the peptide ion\ncharge::UInt8 – Charge of the ion\nmodifier::Float32=PROTON*charge + H2O – Added to the mass of the ion\nisotope::UInt8=UInt8(0) – Diference in the number of isotopes from the monoisotopic ion. \n\nOutput\n\nA Float32 representing the mass-to-charge ratio (m/z) of an ion\n\nNotes\n\nThe modifier argument ought to depend on the kind of ion. For B ions PROTONcharge is appropriate, but for 'y' or precursor ions, PROTONcharge + H2O would be appropriate.\n\nAlgorithm\n\nSum the amino acid residue masses, add modifier + isotope*NEUTRON and then divide the total by the charge. \n\nExamples\n\n#Gets the b6+1 ion MZ\n\njulia> getIonMZ(getResidues(\"PEPTIDE\")[1:6], UInt8(1), modifier = PROTON)\n653.314f0\n\n#Gets the y6+1 ion MZ\n\njulia> getIonMZ(reverse(getResidues(\"PEPTIDE\"))[1:6], UInt8(1))\n703.3144f0\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getIonSeries-Union{Tuple{T}, Tuple{Array{Residue{T}, 1}, UInt8}} where T<:AbstractFloat","page":"Index","title":"Titus.getIonSeries","text":"getIonSeries(residues::Vector{Residue}, charge::UInt8; start::Int = 3, modifier::Float32 = PROTON*charge + H2O, isotope::UInt8 = UInt8(0))\n\nGets the m/z's for an ion series as a Vector{Float32}. \n\nInput\n\nresidues::Vector{Residue}: – List of amino acid residues in the peptide ion\ncharge::UInt8 – Charge of the fragment ions\nstart::Int=3  – Index of first ion the the series to compute.\nisotope::UInt8=UInt8(0) – Diference in the number of isotopes from the monoisotopic ion. \n\nOutput\n\nA Vector{Float32} wich each m/z in the ion series\n\nNotes\n\nThe modifier argument ought to depend on the kind of ion. For a 'b' ion series PROTON*charge is appropriate,\n\nbut for a 'y' ion series, PROTON*charge + H2O would be appropriate. \n\nWill not allow the index of the ion to be equal to or less than the charge. For example,\n\nb2+2 ions could only be calculated in error and are therefore excluded even if start  is set to 2. \n\nIf start exceeds length(residues)-1, then only the N-1 ion is calculated, that is,\n\nthe highest mass ion in the series. \n\nExamples\n\n#Gets the y3+2 through y6+2 ions\n\njulia> getIonSeries(reverse(getResidues(\"PEPTIDE\")), UInt8(2), start = 3)\n4-element view(::Vector{Float32}, 3:6) with eltype Float32:\n 188.58934\n 239.11317\n 287.63956\n 352.16086\n\n#Gets the b4+2 through b6+2 ions\n\njulia> getIonSeries(getResidues(\"PEPTIDE\"), UInt8(2), start = 4, modifier = PROTON*UInt8(2))\n3-element view(::Vector{Float32}, 4:6) with eltype Float32:\n 213.10518\n 269.6472\n 327.16064\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getNearest-Union{Tuple{U}, Tuple{T}, Tuple{Transition, Array{Union{Missing, T}, 1}, Int64}} where {T, U<:AbstractFloat}","page":"Index","title":"Titus.getNearest","text":"getNearest(transition::Transition, masses::Vector{Union{Missing, Float32}}, peak::Int; δ = 0.01)\n\nFinds the peak (index of masses) nearest in mass to the transition but still within the tolerance (getLow(transition)<masses[peak]<getHigh(transition)).  Starts searching at initially supplied peak and increases the index until outside the tolerance. There could be multiple peaks within the tolerance,  and this function selects the one with the lowest mass error to the fragment ion. \n\nInput\n\ntransition::Transition: – Represents a fragment ion\nmasses::Vector{Union{Missing, Float32}} – Mass list from a centroided mass spectrum. MUST BE SORTED IN ASCENDING ORDER. \npeak::Int – An index for a mass in masses\nδ – A mass offset that can be applied to each mass in masses \n\nOutput\n\nAn Int representing the index of the m/z in masses nearest in m/z to that of the fragment m/z. \n\nNotes\n\nIt is assumed that masses is sorted in ascending order. \nIn practice, when called from matchPeaks!, masses[peak] will already be within the fragment m/z tolerance.\n\nUsually there will not be another peak in masses where this is true, but it is possible for multiple peaks to  fall within the tolerance. The purpose of this function is to select the best peak (closes in mass) when this happens.  \n\nAlgorithm\n\nExamples\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getPrecursors-Union{Tuple{Array{Residue{T}, 1}}, Tuple{T}} where T<:AbstractFloat","page":"Index","title":"Titus.getPrecursors","text":"getPrecursors(residues::Vector{Residue}; charges::Vector{UInt8} = UInt8[1,2],isotopes::Vector{UInt8}=UInt8[0],pep_id::UInt32=UInt32(0)) \n\nAlternate constructor for the Precursor struct that Can accept a string representation of a peptide\n\nInput\n\nOutput\n\nNotes\n\n(link to getResidues())\n\n\n\n\n\n","category":"method"},{"location":"#Titus.getYIonModifier-Union{Tuple{UInt8}, Tuple{T}, Tuple{UInt8, T}, Tuple{UInt8, T, T}} where T<:AbstractFloat","page":"Index","title":"Titus.getYIonModifier","text":"getYIonModifier(charge::UInt8)\n\nMass modification that needs to be added to y-ions\n\n\n\n\n\n","category":"method"},{"location":"#Titus.matchPeaks!-Union{Tuple{U}, Tuple{T}, Tuple{Array{FragmentMatch{T}, 1}, Vector{Transition}, Array{Union{Missing, T}, 1}, Array{Union{Missing, T}, 1}, U, UInt32, UInt32}} where {T, U<:AbstractFloat}","page":"Index","title":"Titus.matchPeaks!","text":"function matchPeaks!(matches::Vector{FragmentMatch}, Transitions::Vector{Transition}, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}}, δ::Float64)\n\nFinds the best matching peak in a mass spectrum for each transition/fragment ion supplied if the match is within the fragment tolerance.      Adds each FragmentMatch to matches if not already present. Otherwise,      modifies an existing match (see setFragmentMatch!).\n\nInput\n\nmatches::Vector{FragmentMatch}: – A list representing fragment ions that match peaks in the mass spectrum (masses)\nTransitions::Vector{Transition – A list of fragment ions to search for in the spectrum (masses). MUST BE SORTED IN ASCENDING ORDER BY getMZ(transition)\nmasses::Vector{Union{Missing, Float32}} – Mass list from a centroided mass spectrum. MUST BE SORTED IN ASCENDING ORDER.\nintensities::Vector{Union{Missing, Float32}} – The intensity list from a centroided mass spectrum. Must be the same length as masses\nδ::Float64 – A mass offset that can be applied to each mass in masses\n\nOutput\n\nModifies matches[match] if match is <= lenth(matches). Otherwise adds a new FragmentMatch at matches[match]\n\nNotes\n\nUpdating a match in matches could be useful if researching the same spectra many times at different mass offsets with matchPeaks!    This could be done to calculate a cross correlatoin score for example. If a spectrum is only searched once, then matches should   only be added to matches and existing ones never modified. \nThe fragment tolerance is specified by each Transition. A Transition<:Ion has a field mz::MzFeature which specifies the monoisotopic mass   and also upper and lower bounds (the tolerance). See getLow(ion::Ion) and getHigh(ion::Ion). This is why the user need not supply a fragment tolerance to matchPeaks! \nmasses and intensities contain type unions Union{Missing, Float32}. This method does nothing to check for Missing values, and indeed,   it is assumed that there are none, and the presence of any Missing values will cause an error. The reason for the type union is an idiosyncracy   of the Arrow.jl package and how it implements nested data types in Arrow files. \n\nAlgorithm\n\nGiven a list of fragment ions and a centroided mass spectrum both sorted by m/z, it is efficient to search the spetrum for matches in a \"single pass\" through the spectrum. If there are T transitions and P peaks should be O(T+P). If there are multiple peaks within the tolerance for a given  fragment ion, the peak closest in m/z to the fragment ion is chosen. It is possible to assign the same peak to multiple fragment ions, but  each fragment ion is only assigned to 0 or 1 peaks. \n\nExamples\n\n\n\n\n\n","category":"method"},{"location":"#Titus.matchPeaks-Union{Tuple{U}, Tuple{T}, Tuple{Vector{Transition}, Array{Union{Missing, T}, 1}, Array{Union{Missing, T}, 1}}} where {T, U<:AbstractFloat}","page":"Index","title":"Titus.matchPeaks","text":"function matchPeaks!(matches::Vector{FragmentMatch}, Transitions::Vector{Transition}, masses::Vector{Union{Missing, Float32}}, intensities::Vector{Union{Missing, Float32}}, δ::Float64)\n\nA wrapper for calling matchPeaks at different to search spectra at a list of mass offset.      Each all to matchPeaks Finds the best matching peak in a mass spectrum for each transition/fragment ion supplied if the match is within the fragment tolerance.      Adds each FragmentMatch to matches if not already present. Otherwise, modifies an existing match (see setFragmentMatch!). (see matchPeaks for additional details)\n\nInput\n\nmatches::Vector{FragmentMatch}: – A list representing fragment ions that match peaks in the mass spectrum (masses)\nTransitions::Vector{Transition – A list of fragment ions to search for in the spectrum (masses). MUST BE SORTED IN ASCENDING ORDER BY getMZ(transition)\nmasses::Vector{Union{Missing, Float32}} – Mass list from a centroided mass spectrum. MUST BE SORTED IN ASCENDING ORDER.\nintensities::Vector{Union{Missing, Float32}} – The intensity list from a centroided mass spectrum. Must be the same length as masses\nδ::Float64 – A mass offset that can be applied to each mass in masses\n\nOutput\n\nModifies matches[match] if match is <= lenth(matches). Otherwise adds a new FragmentMatch at matches[match]\n\nNotes\n\nSearching a mass spectrum many times at different mass offsets could be useful for caculating cross correlation scores. \n\nAlgorithm\n\nSee `matchPeaks`\n\nExamples\n\n\n\n\n\n","category":"method"},{"location":"#Titus.setFragmentMatch!-Union{Tuple{T}, Tuple{Array{FragmentMatch{T}, 1}, Int64, Transition, T, T, Int64, UInt32, UInt32}} where T<:AbstractFloat","page":"Index","title":"Titus.setFragmentMatch!","text":"setFragmentMatch!(matches::Vector{FragmentMatch}, match::Int, transition::Transition, mass::Float32, intensity::Float32, peak_ind::Int64)\n\nAdds a FragmentMatch to matches if match is not an index in matches, otherwise, updates the match.\n\nInput\n\nhits::Vector{FragmentMatch}: – Represents a fragment ion\nmatch::Int – Index of the match. Must be <=N+1 where N is length(hits) \nmass::Float32 – m/z of the emperical peak matched to the transition\nintensity::Float32 – intensity of the emperical peak matched to the transition\npeak_ind – unique index of the emperical peak matched to the transition\n\nOutput\n\nModifies matches[match] if match is <= lenth(matches). Otherwise adds a new FragmentMatch at matches[match]\n\nNotes\n\nUpdating a match in matches could be useful if researching the same spectra many times at different mass offsets with matchPeaks!    This could be done to calculate a cross correlatoin score for example. If a spectrum is only searched once, then matches should   only be added to matches and existing ones never modified. \nRecording the peak_ind could be useful, for example, \"chimeric\" scoring of a spectrum from a Vector{FragmentMatch} type. The best scoring precursor would   have fragments matching to known peak_ind. The Vector{FragmentMatch} could be rescored but excluding FragmentMatches corresponding   to those peak_ind's. This would enable a simple chimeric spectra scoring that would not involve completely researching the spectrum (making an additional call to matchPeaks). \n\nAlgorithm\n\nExamples\n\n\n\n\n\n","category":"method"}]
}
