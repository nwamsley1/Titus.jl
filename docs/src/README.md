
## Aims
  Titus aims to simplify preparation and analysis of IS-PRM experiments [1]. Supports analysis of survey runs to characterize internal standard peptides, preparation of methods, peak area ratio estimation and protein-level quantitation. Generates high-quality, informative chromatogram and spectra plots in multi-page pdf's, which are ctrl+f searchable in any standard pdf viewer. Titus is fast, supports multi-threading, and enables analysis of >100 experiments possible in ~1 min including compilation time. Chromatogram plot generation make take a bit longer. Still a work-in-progress and feedback welcome. 
  
## Features
- Cross platform (tested on MacOS and ubuntu)
- Accepts raw ms data in the Apache Arrow format. See the following for cross-platform conversion of Thermo .raw files to the arrow format https://github.com/nwamsley1/ThermoRawFileToParquetConverter.
- Analysis of survey methods. Given a list table of protein-peptide pairs, identifies the best charge state for each precursor (by XTandem hyperscore), the best transitions, and the MS1 peak height. If the survey analyses are split accross multiple experiments, these can be analyzed at once and combined. In addition, can run survey analyses at multiple collision energies/FAIMS CV's to identify the optimum for each analyte. Output is given in a format freindly to XCalibur method editor for Thermo Tribrid instruments.
- Supports variable and fixed modifications defined by regular expressions and includes examples. 
- Estimates peak area ratios using an MM-Estimator (https://github.com/getzze/RobustModels.jl/blob/main/docs/make.jl). Enables accurate par estimation in the pressence of noisy or interfered transitions. High uncertainty in estimation can be used as grounds for exclusion. 
- Summarizaiton of peptide-level quantitation to protein-level quantitation using the MaxLFQ Algorithm without normalization [2].
- Generates a multi-page pdf for each experiment file including chromatogram plots. 
![alt text](https://github.com/nwamsley1/Titus.jl/blob/main/figures/AADLLVNLDPR.png)

## Future Additions/In-Progress 
- Work in progress 
- More complete documentation. At present only basic usage examples are given, and docstrings for most methods are listed in the CI docs but in no particular order. https://documentation.divio.com/
- Robust benchmarking of PAR estimation using robust linear estimators agasint against previously used methods, such as top-N [3] or exclusion of transitions based on cosine-similarity metrics [4].
- If standard curves are available for absolute quantitation, the ability to convert PAR estimates to abundances
- Better plot annotations
- Currently not available as an Julia package. Usability is akward and could use improvement. Not all parameters can be user-defined yet. 
- Some performance issues in "precursors.jl" that should be reasonably easy to resolve. 
- Future support of standard PRM and DIA using prosit libraries and spectral devonvolution?
- Lots of others. Suggest your own....

## Usage

### Survey Run Analyses
```
julia ./src/Routines/PRM/IS-PRM_SURVEY/routine.jl ./data/test.json ./data/parquet/ ./data/NRF2_SIL.txt
```
|Name                |Default| Short        |Description                    |
 |--------------------|-------|-------------|--------------------|
 |params_json||mandatory|Path to a .json file with the parameters (see Configuration)
 |data_dir||mandatory|"Path to a folder with .arrow MS data tables"
 |precursor_list||mandatory|"Path to a tab delimited table of precursors"
 |--make_plots|true|-p|"Whether to make plots. Defaults to `true`"
 |--print_params|false|-s|"Whether to print the parameters from the json. Defaults to `false`"
 
#### Precursor List Example
```
ABCB6	YYNAESYEVER
ABCB6	IDGQDISQVTQASLR
ABCB6	ALNVLVPIFYR
ABHD4	YVSLPNQNK
ADD2	VNVADEVQR
.
.
.
```
### IS-PRM Analysis
```
julia --threads 24 ./src/Routines/PRM/IS-PRM/routine.jl ./data/IS-PRM_TEST.json ./data/parquet ./data/parquet/transition_list.csv
```
|Name                |Default| Short        |Description                    |
 |--------------------|-------|-------------|--------------------|
 |params_json||mandatory|Path to a .json file with the parameters (see Configuration)
 |data_dir||mandatory|"Path to a folder with .arrow MS data tables"
 |transition_list||mandatory|"Path to a tab delimited table of transitions"
 |--make_plots|true|-p|"Whether to make plots. Defaults to `true`"
 |--print_params|false|-s|"Whether to print the parameters from the json. Defaults to `false`"
 
 #### Transition List Example
```
protein_name,sequence,precursor_charge,precursor_isotope,transition_names
ABCB6,YYNAESYEVER[Harg],2,0,y10+1;y7+1;y8+1;y9+1;y6+1
ABCB6,ALNVLVPIFYR[Harg],2,0,y6+1;y8+1;y7+1;b4+1;y9+1
ABHD4,YVSLPNQNK[Hlys],2,0,b4+1;y7+1;y5+2;y3+1;y5+1
.
.
.
```
## Example Outputs

## Configuration files

### IS-PRM-Survey
```
{
    "right_precursor_tolerance": 0.001,
    "left_precursor_tolerance": 0.001,
    "precursor_rt_tolerance": 0.3,
    "b_ladder_start": 3,
    "y_ladder_start": 3,
    "precursor_charges": [2, 3, 4],
    "precursor_isotopes": [0],
    "transition_charges": [1, 2],
    "transition_isotopes": [0],
    "fragment_match_ppm": 40,
    "minimum_fragment_count": 5,
    "fragments_to_select": 5,
    "precursor_rt_window": 0.3,
    "max_variable_mods": 2,
    "fixed_mods":[
                     ["C","C[Carb]"],
                     ["K$","K[Hlys]"],
                     ["R$","R[Harg]"]
    ],
    "variable_mods":
    [],
    "modification_masses":
        {
        "Carb":57.021464,
        "Harg":10.008269,
        "Hlys":8.014199
        },
    "ms_file_conditions":
        {
            "_35NCE_":"35NCE",
            "_40NCE_":"40NCE",
            "GAPDH":"GAPDH"
        }
}
```

### IS-PRM
```
{
    "right_precursor_tolerance": 0.001,
    "left_precursor_tolerance": 0.001,
    "precursor_rt_tolerance": 0.3,
    "b_ladder_start": 3,
    "y_ladder_start": 3,
    "precursor_charges": [2, 3, 4],
    "precursor_isotopes": [0],
    "transition_charges": [1, 2],
    "transition_isotopes": [0],
    "fragment_match_ppm": 40,
    "minimum_fragment_count": 5,
    "fragments_to_select": 5,
    "precursor_rt_window": 0.3,
    "max_variable_mods": 2,
    "fixed_mods":[
                     ["C","C[Carb]"],
                     ["K$","K[Hlys]"],
                     ["R$","R[Harg]"]
    ],
    "variable_mods":
    [],
    "modification_masses":
        {
        "Carb":57.021464,
        "Harg":10.008269,
        "Hlys":8.014199
        },
    "ms_file_conditions":
        {
            "_35NCE_":"35NCE",
            "_40NCE_":"40NCE",
            "GAPDH":"GAPDH"
        }
}
```
## References
<a id="1">[1]</a> 
Gallien S, Kim SY, Domon B. Large-Scale Targeted Proteomics Using Internal Standard Triggered-Parallel Reaction Monitoring (IS-PRM). Mol Cell Proteomics. 2015 Jun;14(6):1630-44. doi: 10.1074/mcp.O114.043968. Epub 2015 Mar 9. PMID: 25755295; PMCID: PMC4458725.
<br>
<a id="1">[2]</a> 
Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. 2014 Sep;13(9):2513-26. doi: 10.1074/mcp.M113.031591. Epub 2014 Jun 17. PMID: 24942700; PMCID: PMC4159666
<br>
<a id="1">[3]</a> 
Stopfer LE, Flower CT, Gajadhar AS, Patel B, Gallien S, Lopez-Ferrer D, White FM. High-Density, Targeted Monitoring of Tyrosine Phosphorylation Reveals Activated Signaling Networks in Human Tumors. Cancer Res. 2021 May 1;81(9):2495-2509. doi: 10.1158/0008-5472.CAN-20-3804. Epub 2021 Jan 28. PMID: 33509940; PMCID: PMC8137532.
<br>
<a id="1">[4]</a> 
Wamsley et al. Targeted proteomic quantitation of NRF2 signaling and predictive biomarkers in HNSCC. https://doi.org/10.1101/2023.03.13.532474 
