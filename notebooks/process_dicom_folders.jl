### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 510a9e53-5c51-476f-9b61-3ced5b9ef968
using DrWatson

# ╔═╡ 5ebd364e-0c4a-4c42-a21c-1dd01c407152
# ╠═╡ show_logs = false
@quickactivate "dynamic"

# ╔═╡ edecc4a9-d370-4567-8945-7879b2bc3871
using PlutoUI, CairoMakie, DICOM

# ╔═╡ 0e5e5563-960c-4525-a18b-d732d88bda61
using Base.Threads

# ╔═╡ 059eeb6f-09a5-4658-b31d-9cec956b359d
TableOfContents()

# ╔═╡ 9bd11631-2645-479a-b9e0-c98fb06c36d2
md"""
## Load DICOMs and paths
"""

# ╔═╡ 20f60d09-92d7-4d06-a076-444cf298114e
function parse_dicoms_threaded(root_folder)
    # Initialize arrays
    dcms = []
    errors = []

    # Locks for thread safety
    dcms_lock = ReentrantLock()
    errors_lock = ReentrantLock()

    # Loop through each file in each directory
    for (root, directories, files) in walkdir(root_folder)
        @threads for file in files
            try
                dcm_path = joinpath(root, file)
                dcm = dcm_parse(dcm_path)

                lock(dcms_lock)  # Acquire lock before modifying shared array
                try
                    push!(dcms, (dcm, dcm_path))
                finally
                    unlock(dcms_lock)  # Always release the lock
                end
            catch e
                lock(errors_lock)  # Acquire lock before modifying shared array
                try
                    push!(errors, e)
                finally
                    unlock(errors_lock)  # Always release the lock
                end
            end
        end
    end

    return dcms, errors
end

# ╔═╡ 8da9ba2c-fdd3-4e37-9935-11b4bff02c79
path = raw"F:\Dynamic Phantom\backup\DICOM\D202308\DD0414"

# ╔═╡ 034bcbdd-735d-414a-854c-cbd1e8085664
dcms, errors = parse_dicoms_threaded(path)

# ╔═╡ 3b512700-a4ed-4098-b6be-3ba347af3333
md"""
## Organize Into Folders
Use header data from the DICOMs to set up the folders
"""

# ╔═╡ 2170163a-a253-4ca7-a04d-cf2b27f91d5c
function organize_and_move_dicoms(dcms, root_folder)
    dcm_dict = Dict{String, Dict}()
    error_paths = []

    for dcm_tuple in dcms
        dcm, dcm_path = dcm_tuple  # Unpack the tuple into the DICOM object and file path

        patient_id = string(get(dcm.meta, tag"Patient Name", "Unknown"))
        series_number = string(get(dcm.meta, tag"Series Number", "Unknown"))
        series_description = string(get(dcm.meta, tag"Series Description", "Unknown"))

        dcm_dict[patient_id] = get(dcm_dict, patient_id, Dict{String, Dict}())
        dcm_dict[patient_id][series_description] = get(dcm_dict[patient_id], series_description, Dict{String, Array}())
        dcm_dict[patient_id][series_description][series_number] = get(dcm_dict[patient_id][series_description], series_number, [])

        # Store both the DICOM object and the file path in the dictionary
        push!(dcm_dict[patient_id][series_description][series_number], (dcm, dcm_path))
    end

    for (patient_id, patient_dict) in dcm_dict
        for (series_description, series_dict) in patient_dict
            for (series_number, dcm_list) in series_dict
                # Create the folder path based on the keys
                patient_path = joinpath(root_folder, patient_id)
                patient_path = rstrip(patient_path)
                folder_path = joinpath(patient_path, series_description, series_number)
                # Create the folder if it doesn't exist
                mkpath(folder_path)

                for (dcm, dcm_path) in dcm_list
                    # Extract the file name from the original path
                    file_name = basename(dcm_path)

                    # Create the new path where the file will be moved
                    new_path = joinpath(folder_path, file_name)

                    try
                        # Move the file to the new location
                        mv(dcm_path, new_path; force=true)
                    catch e
                        push!(error_paths, dcm_path)
                    end
                end
            end
        end
    end
end

# ╔═╡ f19e32a6-18f4-426a-b5a1-68ca5c6bfa9d
root_folder = raw"F:\Dynamic Phantom"

# ╔═╡ f9f16fb4-47a0-4b2b-ad5d-51f219cdf756
organize_and_move_dicoms(dcms, root_folder)

# ╔═╡ Cell order:
# ╠═510a9e53-5c51-476f-9b61-3ced5b9ef968
# ╠═5ebd364e-0c4a-4c42-a21c-1dd01c407152
# ╠═edecc4a9-d370-4567-8945-7879b2bc3871
# ╠═0e5e5563-960c-4525-a18b-d732d88bda61
# ╠═059eeb6f-09a5-4658-b31d-9cec956b359d
# ╟─9bd11631-2645-479a-b9e0-c98fb06c36d2
# ╠═20f60d09-92d7-4d06-a076-444cf298114e
# ╠═8da9ba2c-fdd3-4e37-9935-11b4bff02c79
# ╠═034bcbdd-735d-414a-854c-cbd1e8085664
# ╟─3b512700-a4ed-4098-b6be-3ba347af3333
# ╠═2170163a-a253-4ca7-a04d-cf2b27f91d5c
# ╠═f19e32a6-18f4-426a-b5a1-68ca5c6bfa9d
# ╠═f9f16fb4-47a0-4b2b-ad5d-51f219cdf756
