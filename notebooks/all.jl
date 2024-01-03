### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ a898d895-8335-4a41-a5af-7fe3469580d6
using DrWatson

# ╔═╡ b39a821c-fee6-4f52-8a8a-622136c88332
# ╠═╡ show_logs = false
@quickactivate "dynamic"

# ╔═╡ 98c5ccaf-a05a-4b45-9491-a04a345f3373
using Revise, PlutoUI, CairoMakie, DICOM

# ╔═╡ f088eed8-dbba-4c41-be2f-a38391de6852
using Base.Threads

# ╔═╡ 78f98113-bf70-45f1-97b7-c8c3847aa232
include(srcdir("dicomutils.jl"));

# ╔═╡ 5a3908b9-7dda-4d37-a9fe-6489b3535490
TableOfContents()

# ╔═╡ 2f4ab29e-99f2-491d-b37f-249378d3b428
md"""
# Calibration
"""

# ╔═╡ b97fce2c-9ac5-4c74-8c82-2be751c39e97
path = "/Volumes/USB DISK/DICOM/D202308"

# ╔═╡ 2de3b597-7f46-4ad2-8cad-94394b657b8c
dcms = []

# ╔═╡ 6ae39d7e-81c4-4851-888a-8432d1f6d8d8
for (root, directories, files) in walkdir(path)
    for file in files
        try
            dcm_path = joinpath(root, file)
            dcm = dcm_parse(dcm_path)
            push!(dcms, dcm)
		catch e
			@info e
        end
    end
end

# ╔═╡ 8d9b4d57-9750-4d1e-a86a-c9d6fd7c339d


# ╔═╡ 5d5c63e7-3756-4a08-bbc3-b709ec53391b


# ╔═╡ 222a35ea-d9a6-421b-a1f1-7c913dda9814
# dir = readdir(path)

# ╔═╡ 60830e8e-d557-4c1a-aa70-dfca0b06c787
# dcms = Vector{DICOM.DICOMData}(undef, length(dir_calibration))

# ╔═╡ d214671e-2f95-42e7-a271-72aced943a40
# @threads for idx in eachindex(dir_calibration)
#     dcm_path = joinpath(calibration_path, dir_calibration[idx])
#     dcms[idx] = dcm_parse(dcm_path)
# end

# ╔═╡ f111d863-8756-4ff7-99a4-b28138f466a1
# begin
#     dcm_dict = Dict{String, Dict}()

#     for dcm in dcms
#         patient_id = try
#             dcm.meta[tag"Patient Name"]
#         catch e
#             if isa(e, KeyError)
#                 continue  # skip this iteration if the key does not exist
#             else
#                 rethrow(e)  # if it's a different exception, throw it again
#             end
#         end

#         protocol_name = try
#             dcm.meta[tag"Protocol Name"]
#         catch e
#             if isa(e, KeyError)
#                 continue  # skip this iteration if the key does not exist
#             else
#                 rethrow(e)  # if it's a different exception, throw it again
#             end
#         end

#         series_number = try
#             dcm.meta[tag"Series Number"]
#         catch e
#             if isa(e, KeyError)
#                 continue  # skip this iteration if the key does not exist
#             else
#                 rethrow(e)  # if it's a different exception, throw it again
#             end
#         end

#         if !haskey(dcm_dict, patient_id)
#             dcm_dict[patient_id] = Dict{String, Dict}()
#         end

#         if !haskey(dcm_dict[patient_id], protocol_name)
#             dcm_dict[patient_id][protocol_name] = Dict{Int64, Array}()  # changed key type to Int64
#         end

#         if haskey(dcm_dict[patient_id][protocol_name], series_number)
#             push!(dcm_dict[patient_id][protocol_name][series_number], dcm)
#         else
#             dcm_dict[patient_id][protocol_name][series_number] = [dcm]
#         end
#     end

# end

# ╔═╡ 41fad38f-52b9-4b76-a897-cba44c2fb579
dcm_dict

# ╔═╡ Cell order:
# ╠═a898d895-8335-4a41-a5af-7fe3469580d6
# ╠═b39a821c-fee6-4f52-8a8a-622136c88332
# ╠═98c5ccaf-a05a-4b45-9491-a04a345f3373
# ╠═5a3908b9-7dda-4d37-a9fe-6489b3535490
# ╟─2f4ab29e-99f2-491d-b37f-249378d3b428
# ╠═f088eed8-dbba-4c41-be2f-a38391de6852
# ╠═78f98113-bf70-45f1-97b7-c8c3847aa232
# ╠═b97fce2c-9ac5-4c74-8c82-2be751c39e97
# ╠═2de3b597-7f46-4ad2-8cad-94394b657b8c
# ╠═6ae39d7e-81c4-4851-888a-8432d1f6d8d8
# ╠═8d9b4d57-9750-4d1e-a86a-c9d6fd7c339d
# ╠═5d5c63e7-3756-4a08-bbc3-b709ec53391b
# ╠═222a35ea-d9a6-421b-a1f1-7c913dda9814
# ╠═60830e8e-d557-4c1a-aa70-dfca0b06c787
# ╠═d214671e-2f95-42e7-a271-72aced943a40
# ╠═f111d863-8756-4ff7-99a4-b28138f466a1
# ╠═41fad38f-52b9-4b76-a897-cba44c2fb579
