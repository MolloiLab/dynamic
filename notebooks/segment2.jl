### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fbd6b22d-e6cc-49ca-af42-505ab7d7cb6f
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(temp = true)

	Pkg.add(["Revise", "PlutoUI", "DICOM", "CairoMakie", "LinearAlgebra", "Images", "StatsBase", "Unitful", "Meshes", "Interpolations", "StaticArrays", "CSV", "DataFrames", "IterTools", "Clustering"])
	Pkg.add(url = "https://github.com/Dale-Black/DICOMUtils.jl")
	Pkg.add(url = "https://github.com/Dale-Black/OrthancTools.jl")
	Pkg.develop(path = "/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/dev/julia/ActiveContours.jl")
	
	using Revise, PlutoUI, DICOM, CairoMakie, LinearAlgebra, Images, StatsBase, Unitful, Meshes, Interpolations, StaticArrays, CSV, DataFrames, IterTools, Clustering
	using DICOMUtils, ActiveContours, OrthancTools
end

# ╔═╡ f60bec39-d76f-4b50-9a57-4b8421f7a273
TableOfContents()

# ╔═╡ 68135cc1-3f03-4e17-be3e-b8aed8066067
md"""
# Docs
"""

# ╔═╡ afabde73-b871-4131-a88e-92115c70bcdb
begin
	insert_radii = [1.2, 3.0, 5.0] ./ 2 # mm
	insert_densities = [0.05, 0.1, 0.25, 0.4] # mg/mm^3
	heart_rates = [0, 60, 90] # bpm
	
	scan_types = ["de_80kv", "de_135kv", "se_80kv", "se_120kv", "se_135kv"]
	slice_thicknesses = [0.5, 1.0, 3.0] # mm
	reconstruction_types = ["fbp", "aidr"]
end;

# ╔═╡ 1725e3f6-f542-4a88-8768-d806c196df62
prod = product(insert_radii, insert_densities, heart_rates, scan_types, slice_thicknesses, reconstruction_types);

# ╔═╡ 459247da-d8b2-41c8-8977-91c882753821
begin
	df = DataFrame(collect(prod))
	rename!(df, [:insert_radii, :insert_densities, :heart_rates, :scan_types, :slice_thickness, :reconstruction_types])
	sort!(df, [:insert_radii, :insert_densities, :heart_rates])
end

# ╔═╡ 775ee73b-1658-4fe4-b837-c2562f22f00c
total_accessions = nrow(df) / 30

# ╔═╡ 7311ec8a-f907-431c-86f5-6d902bf26df6
md"""
# Download From Orthanc
"""

# ╔═╡ ed2f601d-2961-4e52-a167-56344c1b4ccc
md"""
## Get Studies
Insert the IP address associated with the Orthanc server into the input box below and then click "Submit". When the code is finished, you can inspect the files by clicking on the dictionary.
"""

# ╔═╡ c22f86ca-0905-491a-bb78-285ff8aef777
# @bind ip_address confirm(TextField(default="128.200.49.26"))

# ╔═╡ f4638a20-671d-4a9b-b9fd-b745378bb41c
# studies_dict = get_all_studies(ip_address)

# ╔═╡ af1ccb54-6778-47cc-bddb-484a80ca8715
# @bind details confirm(download_info("Accession Number", "Series Number(s)", "Instance Number", "Output Directory"))

# ╔═╡ 64532c2f-d4da-4aac-b884-7726495a8d84
# accession_number, series_num, instance_num, output_dir = details

# ╔═╡ 55819ea5-2a37-4e64-951b-9146acf7e93c
md"""
## Get Series
Insert the accession number into the input box above and click "Submit". When the code is finished, you can inspect the files by clicking on the dictionary.
"""

# ╔═╡ 970356fc-f240-4b4e-85e1-d605608b81b5
# series_dict = get_all_series(studies_dict, accession_number, ip_address)

# ╔═╡ 362a18c3-c3b0-4c63-a122-9de073cc0a80
md"""
## Get Instance(s)
You can insert the series number of interest into the input box above and then click "Submit". When the code is finished, you can inspect the files by clicking on the dictionary.
"""

# ╔═╡ 094f7fd0-ef05-43f4-a020-5fb28d1b18fa
# series_num_vec = parse.(Int, split(series_num, ","))

# ╔═╡ c1904341-c714-4304-ad95-3431114e620b
# begin
# 	instances_dicts = []
# 	for i in series_num_vec
# 		instances_dict = get_all_instances(series_dict, string(i), ip_address)
# 		push!(instances_dicts, instances_dict)
# 	end
# end

# ╔═╡ 31903b88-c09c-454c-aecc-467918c39527
instances_dicts

# ╔═╡ 1fddea9e-c127-4f58-bda4-aa5735d7f9be
md"""
# Download DICOM Instance(s)
Type the folder path above, where you want the DICOM files to be saved (or use a temporary directory via `mktempdir()`) in the code cell below. Then type in the instance number that you want to download and click "Submit".
"""

# ╔═╡ 4a1a2b2f-6870-459c-a4e1-2722cf9be36a
# instance_number = parse(Int64, instance_num)

# ╔═╡ c1495c58-5b8c-424c-974e-8c8fd5a6f6f9
# for i in 1:length(instances_dicts)
# 	global output_path = joinpath(output_dir, string(series_num_vec[i]))
# 	if !isdir(output_path)
# 		mkpath(output_path)
# 	end
# 	download_instances(instances_dicts[i], instance_number, output_path, ip_address)
# end

# ╔═╡ fd9ced9d-5c30-4df7-b67b-0fa3de427c2a
# function download_info(acc, ser, inst, save_folder_path)
	
# 	return PlutoUI.combine() do Child
		
# 		inputs = [
# 			md""" $(acc): $(
# 				Child(TextField(default="2781"))
# 			)""",
# 			md""" $(ser): $(
# 				Child(TextField(default="2"))
# 			)""",
# 			md""" $(inst): $(
# 				Child(TextField(default="1"))
# 			)""",
# 			md""" $(save_folder_path): $(
# 				Child(TextField(default="/Users/daleblack/Desktop/dcms")))
# 			)"""
# 		]
		
# 		md"""
# 		#### Scan Details
# 		Input the relevant DICOM information to download the appropriate scans
# 		$(inputs)
# 		"""
# 	end
# end

# ╔═╡ efaef2df-d4f1-4a11-a807-66f4db2b82b7
md"""
## Load DICOMs
"""

# ╔═╡ 10c515ea-18fd-41c1-9fbc-e9fcec58ab54
root = "/Volumes/USB DISK/b"

# ╔═╡ bb45d6bd-d481-4691-9f5d-943872abf690
output_path = joinpath(root, "Cardiac 0.5", "14")

# ╔═╡ 8c83ca3f-e673-40c8-a104-b4aff3f9555a
dcms = dcmdir_parse(output_path)

# ╔═╡ 26677e01-68b4-4a6e-95ac-c848e1abc118
dcm_arr = load_dcm_array(dcms);

# ╔═╡ 356f9aa7-6938-4878-9d82-10b5bd7e1593
header = dcms[1].meta;

# ╔═╡ 1c378cc7-bfb9-46c1-869d-32346bf71ff2
# pixel_size = get_pixel_size(header)

# ╔═╡ 7f0ea094-e04d-48d9-ae17-a0d88fda915e
pixel_size = [0.5, 0.5, 0.5]

# ╔═╡ 71b2125b-d151-430f-be8a-c47c254cebb2
header[tag"PixelSpacing"]

# ╔═╡ 4e4ae00e-9c94-4180-a1e5-74eb2c705c4d
header

# ╔═╡ 48094066-1785-4be7-adbd-36e22a2e72ab
subheader = header[tag"PerFrameFunctionalGroupsSequence"][1].meta

# ╔═╡ 75baf3f7-3770-49be-a7bd-b426e1e1f50f
subheader[tag"PixelSpacing"]

# ╔═╡ 81659d43-242c-4771-8ada-8ef1cc53014e
md"""
## Visualize
"""

# ╔═╡ 0ae11857-bd44-4c76-82af-efe19997d8f9
@bind z1 PlutoUI.Slider(axes(dcm_arr, 3), default=160, show_value=true)

# ╔═╡ 6dc9da67-d8cd-47b0-aa41-c4f1271c08ab
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, z1], colormap=:grays)

	f
end

# ╔═╡ 5a8e0f77-fda5-46cb-81c5-199b23eb67c2
md"""
# Mask Heart Function
"""

# ╔═╡ 3deded4e-5fe1-4086-9e31-8cf6844ed821
function erode_mask(img, num_erosions)
	new_img = img
	i = 0
	while i < num_erosions
		new_img = erode(new_img)
		i += 1
	end
	return new_img
end

# ╔═╡ d40ca853-0015-404f-88bf-03844c4e7b6d
function create_circle_mask(img, centroids, radius)
    # initialize mask with all zeros
    mask = zeros(size(img))

    # define the center of the circle
    x0, y0 = centroids[1], centroids[2]

    # set all pixels inside the circle to 1 and all pixels outside the circle to 0
    for x in axes(img, 1), y in axes(img, 2)
        if ((x - x0)^2 + (y - y0)^2) <= radius^2
            mask[x, y] = 1
        else
            mask[x, y] = 0
        end
    end
    return Bool.(mask)
end

# ╔═╡ dd8bf7a8-1408-4dc0-beaa-dd8d99afc7f5
function centroids_from_mask(mask)
	cc_labels = label_components(mask)
	largest_connected_component, _ = sort(collect(pairs(countmap(cc_labels[cc_labels .!= 0]))), by=x->x[2], rev=true)
	largest_connected_indices = findall(cc_labels .== largest_connected_component[1])

	new_mask = zeros(size(mask))
	for i in largest_connected_indices
		new_mask[i] = 1
	end
	centroids = Int.(round.(component_centroids(label_components(new_mask))[end]))
end

# ╔═╡ e5ebea06-717a-4bf7-ac08-24c06c362cda
heart_cv = chan_vese(dcm_arr);

# ╔═╡ edfe1387-e68b-4b31-97a9-0c05d1fde5ca
centroids = centroids_from_mask(heart_cv)

# ╔═╡ 7f0e576c-5b03-46f4-a8f3-6907dcda62a1
heart_mask = create_circle_mask(dcm_arr[:, :, 3], centroids, 115);

# ╔═╡ fd11d3b8-06fc-451c-9028-b432d0185437
idxs = getindex.(findall(isone, heart_mask), [1 2]);

# ╔═╡ b581f932-0ef0-4ee7-b4ac-038b51cb87fc
md"""
## Visualize
"""

# ╔═╡ 93551f00-906f-4dbc-964e-56872dd2266a
@bind z2 PlutoUI.Slider(axes(dcm_arr, 3), default=130, show_value=true)

# ╔═╡ ce468a2c-85c4-4a9e-94f2-7afff9f0ce5c
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, z2], colormap=:grays)
	scatter!(idxs; markersize=2, color=(:red, 0.5))

	f
end

# ╔═╡ 08d7f3fc-30ad-49cd-97a2-73857d690dd2
md"""
# Inserts
"""

# ╔═╡ b53ebe70-c416-4d07-9bfc-f97962ae4d85
dcm_heart = dcm_arr .* heart_mask;

# ╔═╡ f2f57a29-84ef-4a56-bbb8-bbae5578cf10
@bind d PlutoUI.Slider(axes(dcm_heart, 3), default=130, show_value=true)

# ╔═╡ f0701bba-825c-43c8-bf8d-abd260c008e2
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_heart[:, :, d], colormap=:grays)
	
	f
end

# ╔═╡ 0ddcb24f-170a-43db-bfcc-b9f1f46d8ac5
md"""
## Endpoints
"""

# ╔═╡ 17682dd2-27ff-4a23-a16e-3adb4f1cb0d7
function find_heart_endpoints(dcm_heart, air_threshold=-2000)
	# Find the indices of the elements in the array that are between 30 and 60
	selected_inds = findall(dcm_heart .<= air_threshold)
	
	# Create boolean array from cartesian indices
	local bool_arr4 = zeros(size(dcm_heart))
	for i in selected_inds
		bool_arr4[i] = 1
	end

	local cc_labels 
	while true
		# Use connected component labeling to identify and label all connected components
		cc_labels = label_components(bool_arr4)
		if length(unique(cc_labels)) == 3
			break
		end
		
		bool_arr4 = erode(bool_arr4)
		if (length(unique(cc_labels)) <= 2)
			@warn "Erosion failed; only background is left"
			break
		end
	end
	
	begining_slice = findlast(isone, cc_labels)[3] + 1
	end_slice = findfirst(x->x==2, cc_labels)[3] - 1
	return begining_slice, end_slice
end

# ╔═╡ 53641598-b20a-4186-9251-f63b549f5518
begining_slice, end_slice = find_heart_endpoints(dcm_heart)

# ╔═╡ ec8b5bda-dda5-433f-a686-8128c00e56a1
md"""
## Plane Fitting

**TODO**

Adjust this portion so that the portion of the dcm_heart that is looked at is an eroded version of the heart. This is will ensure that only the air from the insert plane is included and not the outer ring around the heart insert
"""

# ╔═╡ acf73772-3249-49f2-9b41-a2e918261195
function find_heart_plane(dcm_heart, endpoints; air_threshold = -300, std_threshold = 0.25)
	# Find the indices of the elements in the array that are less than
	# air_threshold and filter to exclude beginning and end slices
	remove = [collect(1:endpoints[1])..., collect(endpoints[2]:320)...]
	selected_inds = findall(dcm_heart .<= air_threshold)
	for r in remove
		selected_inds = [i for i in selected_inds if i.I[3] != r]
	end

	selected_inds = getindex.(selected_inds, [1 2 3])

	# Clean points
	centroid = mean(selected_inds, dims=1)
	
	# Calculate the Euclidean distance of each point to the centroid
	distances_to_centroid = [norm(selected_inds[i, :] .- centroid) for i in axes(selected_inds, 1)]
	
	# Calculate mean and standard deviation of the distances
	mean_distance = mean(distances_to_centroid)
	std_distance = std(distances_to_centroid)
	
	# Define the threshold as mean plus 0.5 standard deviations
	threshold = mean_distance + (std_threshold * std_distance)
	
	# Filter the points based on the threshold
	filtered_points = selected_inds[distances_to_centroid .<= threshold, :]

	# Create boolean array from cartesian indices
	bool_arr = zeros(size(dcm_heart))
	for i in eachrow(filtered_points)
		bool_arr[i...] = 1
	end

	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(bool_arr)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	return filtered_points, cc_labels
end

# ╔═╡ b3dc99ff-fad1-4c6f-ab10-c45522adfc21
idx_plane, cc_labels = find_heart_plane(dcm_heart, (begining_slice, end_slice));

# ╔═╡ c7d4a178-6baf-4cd4-b062-42cb17fd3289
function extract_plane_points(cc_labels)
	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# # Find the value with the most occurrences
	most_common_value, _ = sort(collect(pairs(counts)), by=x->x[2], rev=true)

	# # Find the indices of the most common value in the original array
	most_common_indices = findall(cc_labels .== most_common_value[1])
	most_common_indices = getindex.(most_common_indices, [1 2 3])

	pts = most_common_indices[1+100, :], most_common_indices[Int(round(end/2)), :], most_common_indices[end-100, :];
end;

# ╔═╡ dced142e-521a-4748-b66b-db2f69393a6a
centroids

# ╔═╡ a5c4e654-2901-4866-852a-1b3c24244dbd
pts = extract_plane_points(cc_labels)

# ╔═╡ aad3464c-95cc-4d53-893b-5e5359a9ddcc
md"""
## MPR
"""

# ╔═╡ 6250e897-674f-4f02-9a35-f2ace70164e1
function offset_interpolate(itp, plane, i, j, k)
	itp((plane(i, j) + Meshes.normal(plane)*k).coords...)
end

# ╔═╡ 0451e9e3-f572-4120-a102-a111b3616715
function create_mpr(centroids, points, r1, r2)
	offset = Meshes.Point(centroids) - Meshes.Point(first(points)...)
	plane = Meshes.Plane((Meshes.Point.(points...) .+ Ref(offset))...)
	mpr_itp = linear_interpolation(axes(dcm_heart), dcm_heart; extrapolation_bc = 0.0);
	mpr = [offset_interpolate(mpr_itp, plane, i, j, k) for i in -r1:r1, j in -r1:r1, k in -r2:r2]
end

# ╔═╡ 3da7424f-52d1-4b26-bdff-1bfe5e5b2436
mpr = create_mpr(centroids, pts, 512/2, 320/2);

# ╔═╡ 94d19f68-2f83-4237-bd98-7bc0b3c05420
@bind z4 PlutoUI.Slider(axes(mpr, 3), default=centroids[3], show_value=true)

# ╔═╡ 02aa3ce2-8e2b-4e6a-bbca-854b3e5c27fa
heatmap(mpr[:, :, z4], colormap=:grays)

# ╔═╡ 01a047ca-bcae-41e9-86ba-13baf9635055
let
	f = Figure()
	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(mpr[:, :, 160], colormap=:grays)
	scatter!(centroids[1], centroids[2])
	f
end

# ╔═╡ Cell order:
# ╠═fbd6b22d-e6cc-49ca-af42-505ab7d7cb6f
# ╠═f60bec39-d76f-4b50-9a57-4b8421f7a273
# ╟─68135cc1-3f03-4e17-be3e-b8aed8066067
# ╠═afabde73-b871-4131-a88e-92115c70bcdb
# ╠═1725e3f6-f542-4a88-8768-d806c196df62
# ╠═459247da-d8b2-41c8-8977-91c882753821
# ╠═775ee73b-1658-4fe4-b837-c2562f22f00c
# ╟─7311ec8a-f907-431c-86f5-6d902bf26df6
# ╟─ed2f601d-2961-4e52-a167-56344c1b4ccc
# ╠═c22f86ca-0905-491a-bb78-285ff8aef777
# ╠═f4638a20-671d-4a9b-b9fd-b745378bb41c
# ╠═af1ccb54-6778-47cc-bddb-484a80ca8715
# ╠═64532c2f-d4da-4aac-b884-7726495a8d84
# ╟─55819ea5-2a37-4e64-951b-9146acf7e93c
# ╠═970356fc-f240-4b4e-85e1-d605608b81b5
# ╟─362a18c3-c3b0-4c63-a122-9de073cc0a80
# ╠═094f7fd0-ef05-43f4-a020-5fb28d1b18fa
# ╠═c1904341-c714-4304-ad95-3431114e620b
# ╠═31903b88-c09c-454c-aecc-467918c39527
# ╟─1fddea9e-c127-4f58-bda4-aa5735d7f9be
# ╠═4a1a2b2f-6870-459c-a4e1-2722cf9be36a
# ╠═c1495c58-5b8c-424c-974e-8c8fd5a6f6f9
# ╠═fd9ced9d-5c30-4df7-b67b-0fa3de427c2a
# ╟─efaef2df-d4f1-4a11-a807-66f4db2b82b7
# ╠═10c515ea-18fd-41c1-9fbc-e9fcec58ab54
# ╠═bb45d6bd-d481-4691-9f5d-943872abf690
# ╠═8c83ca3f-e673-40c8-a104-b4aff3f9555a
# ╠═26677e01-68b4-4a6e-95ac-c848e1abc118
# ╠═356f9aa7-6938-4878-9d82-10b5bd7e1593
# ╠═1c378cc7-bfb9-46c1-869d-32346bf71ff2
# ╠═7f0ea094-e04d-48d9-ae17-a0d88fda915e
# ╠═71b2125b-d151-430f-be8a-c47c254cebb2
# ╠═4e4ae00e-9c94-4180-a1e5-74eb2c705c4d
# ╠═48094066-1785-4be7-adbd-36e22a2e72ab
# ╠═75baf3f7-3770-49be-a7bd-b426e1e1f50f
# ╟─81659d43-242c-4771-8ada-8ef1cc53014e
# ╟─0ae11857-bd44-4c76-82af-efe19997d8f9
# ╟─6dc9da67-d8cd-47b0-aa41-c4f1271c08ab
# ╟─5a8e0f77-fda5-46cb-81c5-199b23eb67c2
# ╠═3deded4e-5fe1-4086-9e31-8cf6844ed821
# ╠═d40ca853-0015-404f-88bf-03844c4e7b6d
# ╠═dd8bf7a8-1408-4dc0-beaa-dd8d99afc7f5
# ╠═e5ebea06-717a-4bf7-ac08-24c06c362cda
# ╠═edfe1387-e68b-4b31-97a9-0c05d1fde5ca
# ╠═7f0e576c-5b03-46f4-a8f3-6907dcda62a1
# ╠═fd11d3b8-06fc-451c-9028-b432d0185437
# ╟─b581f932-0ef0-4ee7-b4ac-038b51cb87fc
# ╟─93551f00-906f-4dbc-964e-56872dd2266a
# ╟─ce468a2c-85c4-4a9e-94f2-7afff9f0ce5c
# ╟─08d7f3fc-30ad-49cd-97a2-73857d690dd2
# ╠═b53ebe70-c416-4d07-9bfc-f97962ae4d85
# ╟─f2f57a29-84ef-4a56-bbb8-bbae5578cf10
# ╟─f0701bba-825c-43c8-bf8d-abd260c008e2
# ╟─0ddcb24f-170a-43db-bfcc-b9f1f46d8ac5
# ╠═17682dd2-27ff-4a23-a16e-3adb4f1cb0d7
# ╠═53641598-b20a-4186-9251-f63b549f5518
# ╟─ec8b5bda-dda5-433f-a686-8128c00e56a1
# ╠═acf73772-3249-49f2-9b41-a2e918261195
# ╠═b3dc99ff-fad1-4c6f-ab10-c45522adfc21
# ╠═c7d4a178-6baf-4cd4-b062-42cb17fd3289
# ╠═dced142e-521a-4748-b66b-db2f69393a6a
# ╠═a5c4e654-2901-4866-852a-1b3c24244dbd
# ╟─aad3464c-95cc-4d53-893b-5e5359a9ddcc
# ╠═6250e897-674f-4f02-9a35-f2ace70164e1
# ╠═0451e9e3-f572-4120-a102-a111b3616715
# ╠═3da7424f-52d1-4b26-bdff-1bfe5e5b2436
# ╟─94d19f68-2f83-4237-bd98-7bc0b3c05420
# ╟─02aa3ce2-8e2b-4e6a-bbca-854b3e5c27fa
# ╠═01a047ca-bcae-41e9-86ba-13baf9635055
