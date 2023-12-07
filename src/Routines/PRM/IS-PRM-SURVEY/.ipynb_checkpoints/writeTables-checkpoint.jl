function writeTransitionList(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_charge","precursor_isotope","transition_names"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = push!([row[:protein_names], 
                            row[:sequence],
                            row[:precursor_charge],
                            row[:precursor_isotope]], #replace(row[:sequence], r"\[(.*?)\]" => "")
                            join(row[:transition_names],";"))
            write(io, join(data,",")*"\n")
        end
    end
end

function writeIAPIMethod(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_mz","retention_tome","precursor_intensity","condition","transition_mz"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = append!([row[:protein_names], row[:sequence], row[:precursor_mz], row[:retention_time], row[:ms1_peak_height], row[:condition]], row[:transition_mzs])
            write(io, join(data,",")*"\n")
        end
    end
end

function writePrecursorResults(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_mz","retention_time","precursor_intensity","hyperscore","total_ions","file_name"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = [row[:protein_names], row[:sequence], row[:precursor_mz], row[:retention_time], row[:ms1_peak_height], row[:hyperscore], row[:total_ions], row[:file_name]]
            write(io, join(data,",")*"\n")
        end
    end
end