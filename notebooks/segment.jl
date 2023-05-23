### A Pluto.jl notebook ###
# v0.19.25

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

# ╔═╡ 011b0ae2-48d4-4141-8bdb-c94d87ef0a38
using DrWatson

# ╔═╡ c66fda61-be14-4a08-810a-0d75ba617786
# ╠═╡ show_logs = false
@quickactivate "dynamic"

# ╔═╡ f6750b60-00c9-4f89-bc20-0276a2d6ab96
# ╠═╡ show_logs = false
begin
	using Revise, PlutoUI, DICOM, CairoMakie, LinearAlgebra, Images, StatsBase, Unitful, Meshes, Interpolations, MeshViz
	using DICOMUtils, ActiveContours, OrthancTools
end

# ╔═╡ 7e701979-31f6-45f2-a4d6-6928d93ad042
TableOfContents()

# ╔═╡ 2330d887-782e-4da1-ba3c-f27365ab3203
md"""
# Download From Orthanc
"""

# ╔═╡ ce6fc772-4353-4e34-ad46-c5664c8a431d
md"""
## Get Studies
Insert the IP address associated with the Orthanc server into the input box below and then click "Submit". When the code is finished, you can inspect the files by clicking on the dictionary.
"""

# ╔═╡ 29398382-062e-4904-a09b-487ce65eb3da
@bind ip_address confirm(TextField(default="128.200.49.26"))

# ╔═╡ 7624c25a-7dfa-4805-a014-40a61001c72b
studies_dict = get_all_studies(ip_address)

# ╔═╡ 31701dba-2979-4e86-aebf-d1ebb0306807
md"""
## Get Series
Insert the accession number into the input box above and click "Submit". When the code is finished, you can inspect the files by clicking on the dictionary.
"""

# ╔═╡ 28ce38bf-89e3-40b3-9470-03aeb5bfbfb1
md"""
## Get Instance(s)
You can insert the series number of interest into the input box above and then click "Submit". When the code is finished, you can inspect the files by clicking on the dictionary.
"""

# ╔═╡ a0b0cf10-76cc-4fa0-abc6-2cdbfb880617
md"""
# Download DICOM Instance(s)
Type the folder path above, where you want the DICOM files to be saved (or use a temporary directory via `mktempdir()`) in the code cell below. Then type in the instance number that you want to download and click "Submit".
"""

# ╔═╡ a257fd92-ebe0-4551-a83c-0db89d54dbb4
function download_info(acc, ser, inst, save_folder_path)
	
	return PlutoUI.combine() do Child
		
		inputs = [
			md""" $(acc): $(
				Child(TextField(default="2475"))
			)""",
			md""" $(ser): $(
				Child(TextField(default="2"))
			)""",
			md""" $(inst): $(
				Child(TextField(default="1"))
			)""",
			md""" $(save_folder_path): $(
				Child(TextField(default="/Users/daleblack/Desktop/dcms")))
			)"""
		]
		
		md"""
		#### Scan Details
		Input the relevant DICOM information to download the appropriate scans
		$(inputs)
		"""
	end
end

# ╔═╡ 1ee500b3-7c50-4bd8-bf32-5fee3a373db2
@bind details confirm(download_info("Accession Number", "Series Number(s)", "Instance Number", "Output Directory"))

# ╔═╡ 0ece2fe8-1bd2-45a4-a5dd-94501f553824
accession_number, series_num, instance_num, output_dir = details

# ╔═╡ 98e6b66a-551b-495a-92fe-48d825db7ae5
series_dict = get_all_series(studies_dict, accession_number, ip_address)

# ╔═╡ 6f36409a-3dfc-4a2e-977c-05953ff3a7dc
series_num_vec = parse.(Int, split(series_num, ","))

# ╔═╡ 812a0307-b92f-4c4c-91c8-769a1bb5d6b1
begin
	instances_dicts = []
	for i in series_num_vec
		instances_dict = get_all_instances(series_dict, string(i), ip_address)
		push!(instances_dicts, instances_dict)
	end
end

# ╔═╡ b1f5d8bd-a197-490e-89ab-c3258c4ac2f6
instances_dicts

# ╔═╡ 60de75b1-791c-446b-b99b-beb799a1a9f4
instance_number = parse(Int64, instance_num)

# ╔═╡ 081e437e-dc99-421a-a4b4-735d47cdb77d
for i in 1:length(instances_dicts)
	global output_path = joinpath(output_dir, string(series_num_vec[i]))
	if !isdir(output_path)
		mkpath(output_path)
	end
	download_instances(instances_dicts[i], instance_number, output_path, ip_address)
end

# ╔═╡ 340a1e3d-9d29-4afe-88d3-2e38da413595
md"""
## Load DICOMs
"""

# ╔═╡ b8a51c87-90cf-4259-9061-4e8dc8e0632c
dcms = dcmdir_parse(output_path)

# ╔═╡ b38f92ac-ad86-404b-bf33-118d7f886b36
dcm_arr = load_dcm_array(dcms);

# ╔═╡ 33eae10c-3777-4c53-8245-72932a28069d
header = dcms[1].meta;

# ╔═╡ 3ad9eb7c-9241-4a06-8308-b6ba33aa492f
pixel_size = get_pixel_size(header)

# ╔═╡ 6321ab34-df67-410a-8776-850b4b456058
md"""
## Visualize
"""

# ╔═╡ 3fd05c3f-eb99-4d26-8b55-eb96b1513bc1
@bind a PlutoUI.Slider(axes(dcm_arr, 3), default=160, show_value=true)

# ╔═╡ 496773b5-a3bd-4ed7-aba5-7effe8cd958b
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, a], colormap=:grays)

	f
end

# ╔═╡ 5a76fded-55c3-42bc-bfa5-988e0e5c6f53
md"""
# Mask Heart Function
"""

# ╔═╡ 150f915c-8c03-4fbc-99be-9b88927c60a9
function erode_mask(img, num_erosions)
	new_img = img
	i = 0
	while i < num_erosions
		new_img = erode(new_img)
		i += 1
	end
	return new_img
end

# ╔═╡ 7cf13428-0e4d-48bd-90e7-239a5e1e610c
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

# ╔═╡ b85c9228-1317-484a-abd4-f1ac725352c2
function mask_heart(
	dcm_array, radius1, radius2;
	limits=(40, 70), mu=0.25, lambda1=1, lambda2=1, tol=1e-3, max_iter=200, dt=0.5
)

	## -- Filter Heart Tissue --##
	# Find the indices of the elements in the array that are between 30 and 60
	selected_inds = findall((limits[1] .<= dcm_array) .& (dcm_array.<= limits[2]))

	# Create boolean array from cartesian indices
	bool_arr = zeros(size(dcm_array))
	for i in selected_inds
		bool_arr[i] = 1
	end

	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(bool_arr)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# # Find the value with the most occurrences
	most_common_value, _ = sort(collect(pairs(counts)), by=x->x[2], rev=true)
	
	# # Find the indices of the most common value in the original array
	most_common_indices = findall(cc_labels .== most_common_value[1])

	# Create boolean array from new cartesian indices
	bool_arr2 = zeros(size(dcm_array))
	for i in most_common_indices
		bool_arr2[i] = 1
	end
	
	centroids = Int.(round.(component_centroids(label_components(bool_arr2))[2]))
	dcm_slice = dcm_arr[:, :, centroids[3]]
	circle_mask = create_circle_mask(dcm_slice, centroids[1:2], radius1)

	## -- Chan Vese --##
	checkerboard = init_checkerboard(size(dcm_slice), 5);
	checkerboard = circle_mask .* checkerboard;
	for i in axes(checkerboard, 1)
		for j in axes(checkerboard, 2)
			if checkerboard[i, j] == 0
				checkerboard[i, j] = -1
			end
		end
	end

	segmentation, phi, energies = chan_vese(
		dcm_slice, checkerboard; 
		mu=mu, lambda1=lambda1, lambda2=lambda2, tol=tol, max_iter=max_iter, dt=dt);

	while true
		segmentation = erode(segmentation)
		cc = label_components(segmentation)
		if length(unique(cc)) == 2
			break
		end
		if (length(unique(cc)) <= 1)
			@warn "Erosion failed; only background is left"
			break
		end
	end
	
	segmentation_lbl = label_components(segmentation)
	centroids2 = component_centroids(segmentation_lbl)[2]

	heart_mask = create_circle_mask(dcm_slice, centroids2, radius2)
	return heart_mask, centroids
end

# ╔═╡ b1f2ab09-7fd4-4057-9873-fba70ee4f326
begin
	heart_mask2, centroids = mask_heart(dcm_arr, 130, 70)
	idxs2 = findall(isone, heart_mask2)
	idxs2 = getindex.(idxs2, [1 2])
end;

# ╔═╡ bb4bcf58-7b2e-4719-9d73-f863752c6e59
md"""
## Visualize
"""

# ╔═╡ 430312af-0d9f-4f29-a823-c8230ea6167d
@bind c PlutoUI.Slider(axes(dcm_arr, 3), default=130, show_value=true)

# ╔═╡ 7e3e3519-6851-40c7-b1c4-bfe54f40f355
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, c], colormap=:grays)
	scatter!(idxs2; markersize=2, color=(:red, 0.5))

	f
end

# ╔═╡ 27aa116c-ce51-4b9c-8651-573bd6ff431b
md"""
# Inserts
"""

# ╔═╡ 6cde7a0d-3ec0-4bf2-8e2f-ddb3f877f3a4
dcm_heart = dcm_arr .* heart_mask2;

# ╔═╡ 76c8cdd4-d189-47b1-8780-f5da614f0572
@bind d PlutoUI.Slider(axes(dcm_heart, 3), default=130, show_value=true)

# ╔═╡ a5a087a1-6be0-4a9d-812b-4cacab3717a2
let
	f = Figure()

	CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_heart[:, :, d], colormap=:grays)
	
	f
end

# ╔═╡ 26b87f0e-6588-4dcc-b180-0b8920a6736e
md"""
## Endpoints
"""

# ╔═╡ ae417bdb-34b6-4641-b17e-3cb38b45bff4
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

# ╔═╡ d5fab0cf-1f5d-4651-a535-dcf6beaed0de
begining_slice, end_slice = find_heart_endpoints(dcm_heart)

# ╔═╡ 0272fdba-8045-462f-b087-87fdc5330e6b
md"""
## Plane Fitting
"""

# ╔═╡ 92c189d8-8120-4fec-b627-aa44191dc9f5
function find_heart_plane(dcm_heart, endpoints, air_threshold=-300)
	# Find the indices of the elements in the array that are less than
	# air_threshold and filter to exclude beginning and end slices
	remove = [collect(1:endpoints[1])..., collect(endpoints[2]:320)...]
	selected_inds = findall(dcm_heart .<= air_threshold)
	for r in remove
		selected_inds = [i for i in selected_inds if i.I[3] != r]
	end

	# Create boolean array from cartesian indices
	bool_arr = zeros(size(dcm_heart))
	for i in selected_inds
		bool_arr[i] = 1
	end

	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(bool_arr)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	return selected_inds, cc_labels
end

# ╔═╡ 2fb08833-4fda-4143-aebe-ce9d12e5886c
begin
	idx_plane, cc_labels = find_heart_plane(dcm_heart, (begining_slice, end_slice))
	idx_plane = getindex.(idx_plane, [1 2 3])
end

# ╔═╡ 2aa24477-adc5-45b0-a103-f7fbe577cab0
begin
	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# # Find the value with the most occurrences
	most_common_value, _ = sort(collect(pairs(counts)), by=x->x[2], rev=true)

	# # Find the indices of the most common value in the original array
	most_common_indices = findall(cc_labels .== most_common_value[1])
	most_common_indices = getindex.(most_common_indices, [1 2 3])
end

# ╔═╡ 1a5ed8aa-1f34-4cf6-b1e3-cc1a55fd123d
pts = most_common_indices[1+100, :], most_common_indices[Int(round(end/2)), :], most_common_indices[end-100, :];

# ╔═╡ 510683e8-10fc-4d33-b3c3-a9b7b14e5335
@bind e PlutoUI.Slider(axes(dcm_heart, 3), default=125, show_value=true)

# ╔═╡ 4c971b36-e6ce-4f38-aced-a15d55312a99
let
	idx_plane_2d = findall(idx_plane[:, 3] .== e)
	f, l = first(idx_plane_2d), last(idx_plane_2d)
	idx_plane_2d = idx_plane[f:l, 1:2]
	
	f = Figure()
	
	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(dcm_heart[:, :, e], colormap=:grays)
	scatter!(idx_plane_2d; markersize=2, color=:red)

	f
end

# ╔═╡ 0f11603a-c90f-4130-89fd-9575328871f8
md"""
## MPR
"""

# ╔═╡ ca3ba794-0ae5-43d6-a93a-00f4e625be8f
offset = Meshes.Point(centroids) - Meshes.Point(first(pts));

# ╔═╡ 67735130-1469-4acd-8c44-8b2adf37383c
plane = Meshes.Plane((Meshes.Point.(pts) .+ Ref(offset))...);

# ╔═╡ 8da676da-7175-4a25-a3f1-2fcfa4c60110
offset_interpolate(itp, plane, i, j, k) = itp((plane(i, j) + Meshes.normal(plane)*k).coords...);

# ╔═╡ ec18dafd-fb16-4d84-9f43-dfdb7d7f837f
mpr_itp = linear_interpolation(axes(dcm_heart), dcm_heart; extrapolation_bc = 0.0);

# ╔═╡ f2911ebb-6a75-4f2a-b9c8-e9a8dfb9436a
r, r2 = 512/2, 320/2;

# ╔═╡ 83a525a1-1ec5-4174-adc7-b19a2890c37b
mpr = [offset_interpolate(mpr_itp, plane, i, j, k) for i in -r:r, j in -r:r, k in -r2:r2];

# ╔═╡ 8f6a2277-504d-4009-ac15-cca75b1e6575
@bind f PlutoUI.Slider(axes(mpr, 3), default=centroids[3], show_value=true)

# ╔═╡ b53a9ab0-38e3-4144-990a-55f1f5cdd33d
heatmap(mpr[:, :, f], colormap=:grays)

# ╔═╡ 17277b11-1c96-4ed4-9acb-be25403fe2aa
md"""
## Cavity centerpoints
"""

# ╔═╡ 5bfbe465-2753-4054-b038-eca14e70c089
mpr_slice = mpr[:, :, centroids[3]];

# ╔═╡ 9fb78b86-c9ea-4247-8827-374ecdf1a771
thresh = 200

# ╔═╡ 73c2b825-9de4-4559-98e1-81eba0f57af0
mpr_slice_thresh = mpr_slice .> thresh;

# ╔═╡ db855edd-2f3e-42b4-b8c9-64e949669716
begin
	# Use connected component labeling to identify and label all connected components
	cc_labels2 = label_components(mpr_slice_thresh)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts2 = countmap(cc_labels2[cc_labels2 .!= 0])
	
	# # Find the value with the most occurrences
	most_common_value_a, most_common_value_b, most_common_value_c, _ = sort(collect(pairs(counts2)), by=x->x[2], rev=true)
end

# ╔═╡ b85e535c-d2dd-498a-9bfa-c8ac3965167f
most_common_value_a, most_common_value_b, most_common_value_c

# ╔═╡ a4d0e49a-a962-42d8-a278-32fea583281c
begin
	# # Find the indices of the most common value in the original array
	most_common_indices_a = findall(cc_labels2 .== most_common_value_a[1])

	# Create boolean array from new cartesian indices
	bool_arr_a = zeros(size(mpr_slice))
	for i in most_common_indices_a
		bool_arr_a[i] = 1
	end
	centroids_a = Int.(round.(component_centroids(label_components(bool_arr_a))[2]))
	box_a = component_boxes(label_components(bool_arr_a))

	# Find the indices of the most common value in the original array
	most_common_indices_b = findall(cc_labels2 .== most_common_value_b[1])

	# Create boolean array from new cartesian indices
	bool_arr_b = zeros(size(mpr_slice))
	for i in most_common_indices_b
		bool_arr_b[i] = 1
	end
	centroids_b = Int.(round.(component_centroids(label_components(bool_arr_b))[2]))

	# # Find the indices of the most common value in the original array
	most_common_indices_c = findall(cc_labels2 .== most_common_value_c[1])

	# Create boolean array from new cartesian indices
	bool_arr_c= zeros(size(mpr_slice))
	for i in most_common_indices_c
		bool_arr_c[i] = 1
	end
	centroids_c = Int.(round.(component_centroids(label_components(bool_arr_c))[2]))
end

# ╔═╡ 4121b935-7b01-4176-9cb9-115e61e401e5
heatmap(label_components(bool_arr_a))

# ╔═╡ c0b40d91-d9fa-44e6-8f90-26a2a14a0bab
box_a[2]

# ╔═╡ 18430fe7-5612-43a0-8ae1-fb240f97c269
box_a[2][1]

# ╔═╡ 565ebd37-cd13-49ad-88ef-e3b9d0de3bec
centers_a, centers_b, centers_c = [centroids_a..., centroids[3]], [centroids_b..., centroids[3]], [centroids_c..., centroids[3]]

# ╔═╡ 35f9fea2-d84e-496d-b433-cd4aaebac40c
pt_a, pt_b, pt_c = Meshes.Point(centers_a), Meshes.Point(centers_b), Meshes.Point(centers_c)

# ╔═╡ 6483fd22-acdb-4325-8d67-bd54c08c9a27
let
	msize = 5
	f = Figure()

	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(mpr_slice; colormap=:grays)
	scatter!(box_a[2][1], markersize=msize, color=:red)
	scatter!(box_a[2][2], markersize=msize, color=:red)
	scatter!(centroids_a, markersize=msize, color=:purple)
	scatter!(centroids_b, markersize=msize, color=:blue)
	scatter!(centroids_c, markersize=msize, color=:green)

	f
end

# ╔═╡ 200c67b9-bff2-456c-9b74-0ca1c31ee964
md"""
## Cavity cylinders
"""

# ╔═╡ 17e65e3a-5a10-4004-ab11-4018167db3e6
begin
	vec_a = pt_a - pt_c
	plane_a = Meshes.Plane(pt_a, vec_a)

	vec_c = pt_c - pt_a
	plane_c = Meshes.Plane(pt_c, vec_c)
end

# ╔═╡ f7e370da-7542-4277-84a5-a9c64efaf27b
cyl = Cylinder(10.0, plane_a, plane_c)

# ╔═╡ 51e62cb9-a957-4c1b-83d8-7ed7b0cc87b9
viz(cyl)

# ╔═╡ Cell order:
# ╠═011b0ae2-48d4-4141-8bdb-c94d87ef0a38
# ╠═c66fda61-be14-4a08-810a-0d75ba617786
# ╠═f6750b60-00c9-4f89-bc20-0276a2d6ab96
# ╠═7e701979-31f6-45f2-a4d6-6928d93ad042
# ╟─2330d887-782e-4da1-ba3c-f27365ab3203
# ╟─ce6fc772-4353-4e34-ad46-c5664c8a431d
# ╟─29398382-062e-4904-a09b-487ce65eb3da
# ╠═7624c25a-7dfa-4805-a014-40a61001c72b
# ╟─1ee500b3-7c50-4bd8-bf32-5fee3a373db2
# ╠═0ece2fe8-1bd2-45a4-a5dd-94501f553824
# ╟─31701dba-2979-4e86-aebf-d1ebb0306807
# ╠═98e6b66a-551b-495a-92fe-48d825db7ae5
# ╟─28ce38bf-89e3-40b3-9470-03aeb5bfbfb1
# ╠═6f36409a-3dfc-4a2e-977c-05953ff3a7dc
# ╠═812a0307-b92f-4c4c-91c8-769a1bb5d6b1
# ╠═b1f5d8bd-a197-490e-89ab-c3258c4ac2f6
# ╟─a0b0cf10-76cc-4fa0-abc6-2cdbfb880617
# ╠═60de75b1-791c-446b-b99b-beb799a1a9f4
# ╠═081e437e-dc99-421a-a4b4-735d47cdb77d
# ╟─a257fd92-ebe0-4551-a83c-0db89d54dbb4
# ╟─340a1e3d-9d29-4afe-88d3-2e38da413595
# ╠═b8a51c87-90cf-4259-9061-4e8dc8e0632c
# ╠═b38f92ac-ad86-404b-bf33-118d7f886b36
# ╠═33eae10c-3777-4c53-8245-72932a28069d
# ╠═3ad9eb7c-9241-4a06-8308-b6ba33aa492f
# ╟─6321ab34-df67-410a-8776-850b4b456058
# ╟─3fd05c3f-eb99-4d26-8b55-eb96b1513bc1
# ╟─496773b5-a3bd-4ed7-aba5-7effe8cd958b
# ╟─5a76fded-55c3-42bc-bfa5-988e0e5c6f53
# ╠═150f915c-8c03-4fbc-99be-9b88927c60a9
# ╠═7cf13428-0e4d-48bd-90e7-239a5e1e610c
# ╠═b85c9228-1317-484a-abd4-f1ac725352c2
# ╠═b1f2ab09-7fd4-4057-9873-fba70ee4f326
# ╟─bb4bcf58-7b2e-4719-9d73-f863752c6e59
# ╟─430312af-0d9f-4f29-a823-c8230ea6167d
# ╟─7e3e3519-6851-40c7-b1c4-bfe54f40f355
# ╟─27aa116c-ce51-4b9c-8651-573bd6ff431b
# ╠═6cde7a0d-3ec0-4bf2-8e2f-ddb3f877f3a4
# ╟─76c8cdd4-d189-47b1-8780-f5da614f0572
# ╟─a5a087a1-6be0-4a9d-812b-4cacab3717a2
# ╟─26b87f0e-6588-4dcc-b180-0b8920a6736e
# ╠═ae417bdb-34b6-4641-b17e-3cb38b45bff4
# ╠═d5fab0cf-1f5d-4651-a535-dcf6beaed0de
# ╟─0272fdba-8045-462f-b087-87fdc5330e6b
# ╠═92c189d8-8120-4fec-b627-aa44191dc9f5
# ╠═2fb08833-4fda-4143-aebe-ce9d12e5886c
# ╠═2aa24477-adc5-45b0-a103-f7fbe577cab0
# ╠═1a5ed8aa-1f34-4cf6-b1e3-cc1a55fd123d
# ╟─510683e8-10fc-4d33-b3c3-a9b7b14e5335
# ╟─4c971b36-e6ce-4f38-aced-a15d55312a99
# ╟─0f11603a-c90f-4130-89fd-9575328871f8
# ╠═ca3ba794-0ae5-43d6-a93a-00f4e625be8f
# ╠═67735130-1469-4acd-8c44-8b2adf37383c
# ╠═8da676da-7175-4a25-a3f1-2fcfa4c60110
# ╠═ec18dafd-fb16-4d84-9f43-dfdb7d7f837f
# ╠═f2911ebb-6a75-4f2a-b9c8-e9a8dfb9436a
# ╠═83a525a1-1ec5-4174-adc7-b19a2890c37b
# ╟─8f6a2277-504d-4009-ac15-cca75b1e6575
# ╟─b53a9ab0-38e3-4144-990a-55f1f5cdd33d
# ╟─17277b11-1c96-4ed4-9acb-be25403fe2aa
# ╠═5bfbe465-2753-4054-b038-eca14e70c089
# ╠═9fb78b86-c9ea-4247-8827-374ecdf1a771
# ╠═73c2b825-9de4-4559-98e1-81eba0f57af0
# ╠═db855edd-2f3e-42b4-b8c9-64e949669716
# ╠═b85e535c-d2dd-498a-9bfa-c8ac3965167f
# ╠═a4d0e49a-a962-42d8-a278-32fea583281c
# ╠═4121b935-7b01-4176-9cb9-115e61e401e5
# ╠═c0b40d91-d9fa-44e6-8f90-26a2a14a0bab
# ╠═18430fe7-5612-43a0-8ae1-fb240f97c269
# ╠═565ebd37-cd13-49ad-88ef-e3b9d0de3bec
# ╠═35f9fea2-d84e-496d-b433-cd4aaebac40c
# ╟─6483fd22-acdb-4325-8d67-bd54c08c9a27
# ╟─200c67b9-bff2-456c-9b74-0ca1c31ee964
# ╠═17e65e3a-5a10-4004-ab11-4018167db3e6
# ╠═f7e370da-7542-4277-84a5-a9c64efaf27b
# ╠═51e62cb9-a957-4c1b-83d8-7ed7b0cc87b9
