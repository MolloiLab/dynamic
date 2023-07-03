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

# ╔═╡ 011b0ae2-48d4-4141-8bdb-c94d87ef0a38
using DrWatson

# ╔═╡ c66fda61-be14-4a08-810a-0d75ba617786
# ╠═╡ show_logs = false
@quickactivate "dynamic"

# ╔═╡ f6750b60-00c9-4f89-bc20-0276a2d6ab96
# ╠═╡ show_logs = false
begin
	using Revise, PlutoUI, DICOM, CairoMakie, LinearAlgebra, Images, StatsBase, Unitful, Meshes, Interpolations, StaticArrays
	using DICOMUtils, ActiveContours, OrthancTools
end

# ╔═╡ b6f3279a-1c31-4c42-bac2-1adc9ec317b4
begin
	using CalciumScoring
	using CalciumScoring: score
end

# ╔═╡ 2d0f900a-85e2-4956-b08e-9f85e1d68a60
include(srcdir("masks.jl"));

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
				Child(TextField(default="2751"))
			)""",
			md""" $(ser): $(
				Child(TextField(default="5"))
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

# ╔═╡ 4edf70c6-a181-494d-b996-8bb8a44f1f32
# output_path = "/Users/daleblack/Desktop/dcms/5"

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
	heart_mask2, centroids = mask_heart(dcm_arr, 130, 90)
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

# ╔═╡ 4b50a86e-fb36-4d26-9460-111f0465dce7
function offset_interpolate(itp, plane, i, j, k)
	itp((plane(i, j) + Meshes.normal(plane)*k).coords...)
end

# ╔═╡ 6edfa2b7-63ae-4315-ac9f-caf96f0f51ed
Meshes.Point(first(pts)...)

# ╔═╡ 32d8eec0-892e-455f-9815-26eccc4dd3d6
function create_mpr(centroids, points, r1, r2)
	offset = Meshes.Point(centroids) - Meshes.Point(first(points))
	plane = Meshes.Plane((Meshes.Point.(points) .+ Ref(offset))...)
	mpr_itp = linear_interpolation(axes(dcm_heart), dcm_heart; extrapolation_bc = 0.0);
	mpr = [offset_interpolate(mpr_itp, plane, i, j, k) for i in -r1:r1, j in -r1:r1, k in -r2:r2]
end

# ╔═╡ d3282a5e-1a7a-4537-b7f1-cc4a7c99c464
mpr = create_mpr(centroids, pts, 512/2, 320/2);

# ╔═╡ 24c102fd-1fc0-4c31-9ad1-d9af691159b0
centroids

# ╔═╡ 8f6a2277-504d-4009-ac15-cca75b1e6575
@bind f PlutoUI.Slider(axes(mpr, 3), default=centroids[3], show_value=true)

# ╔═╡ b53a9ab0-38e3-4144-990a-55f1f5cdd33d
heatmap(mpr[:, :, f], colormap=:grays)

# ╔═╡ aa8d9f5a-c272-4dba-894c-653b92501bb3
let
	f = Figure()
	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(mpr[:, :, 158], colormap=:grays)
	scatter!(centroids[1], centroids[2])
	f
end

# ╔═╡ 200c67b9-bff2-456c-9b74-0ca1c31ee964
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 97739511-af5a-46f9-a251-09c1e13ce758
function get_insert_centers(mpr, threshold)
	z = div(size(mpr, 3), 2)
	mpr_slice = mpr[:, :, z]
	mpr_slice_thresh = mpr_slice .> threshold
	
	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(mpr_slice_thresh)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# Find the value with the most occurrences
	most_common_value_a, most_common_value_b = sort(collect(pairs(counts)), by=x->x[2], rev=true)

	# Find the indices of the most common value in the original array
	most_common_indices_a = findall(cc_labels .== most_common_value_a[1])

	# Create boolean array from new cartesian indices
	bool_arr_a = zeros(size(mpr_slice))
	for i in most_common_indices_a
		bool_arr_a[i] = 1
	end
	centroids_a = Int.(round.(component_centroids(label_components(bool_arr_a))[2]))
	box_a = component_boxes(label_components(bool_arr_a))

	# Find the indices of the most common value in the original array
	most_common_indices_b = findall(cc_labels .== most_common_value_b[1])

	# Create boolean array from new cartesian indices
	bool_arr_b = zeros(size(mpr_slice))
	for i in most_common_indices_b
		bool_arr_b[i] = 1
	end
	centroids_b = Int.(round.(component_centroids(label_components(bool_arr_b))[2]))

	centers_a, centers_b = [centroids_a..., z], [centroids_b..., z]
	return centers_a, centers_b
	
end

# ╔═╡ 216c4a93-4acc-4c3e-9856-40ed1daf4611
centers_a, centers_b = get_insert_centers(mpr, 200);

# ╔═╡ 6483fd22-acdb-4325-8d67-bd54c08c9a27
let
	msize = 5
	f = Figure()

	ax = CairoMakie.Axis(f[1, 1])
	z = div(size(mpr, 3), 2)
	heatmap!(mpr[:, :, z]; colormap=:grays)
	scatter!(centers_a[1], centers_a[2], markersize=msize, color=:purple)
	scatter!(centers_b[1], centers_b[2], markersize=msize, color=:blue)

	f
end

# ╔═╡ 649f943d-8838-4258-a108-4de678561d83
# Modify the in_cylinder function to accept Static Vectors
function _in_cylinder(pt::SVector{3, Int}, pt1::SVector{3, Float64}, pt2::SVector{3, Float64}, radius)
    v = pt2 - pt1
    w = pt - pt1

    # Compute the dot product
    c1 = dot(w, v)
    if c1 <= 0
        return norm(w) <= radius
    end

    c2 = dot(v, v)
    if c2 <= c1
        return norm(pt - pt2) <= radius
    end

    # Compute the perpendicular distance
    b = c1 / c2
    pb = pt1 + b * v
    return norm(pt - pb) <= radius
end

# ╔═╡ 2e72c9bb-e41c-493f-b62f-9cc60531382a
function create_cylinder(array, pt1, pt2, radius, offset)
    # Convert the points to static arrays
    pt1 = SVector{3, Float64}(pt1)
    pt2 = SVector{3, Float64}(pt2)

    # Compute the unit vector in the direction from pt1 to pt2
    direction = normalize(pt2 - pt1)

    # Adjust the endpoints of the cylinder by the offset
    pt1 = pt1 - offset * direction
    pt2 = pt2 + offset * direction

    # Initialize the 3D array
    cylinder = zeros(Int, size(array)...)
    # Iterate over the 3D array
    for k in axes(cylinder, 3)
        for j in axes(cylinder, 2)
            for i in axes(cylinder, 1)
                # Create a static vector for the current point
                pt = SVector{3, Int}(i, j, k)

                # Check if the current point is inside the cylinder
                if _in_cylinder(pt, pt1, pt2, radius)
                    cylinder[i, j, k] = 1
                end
            end
        end
    end
    return Bool.(cylinder)
end

# ╔═╡ 6b1ee364-7611-4fde-9354-aeee659a127f
cylinder = create_cylinder(mpr, centers_a, centers_b, 8, -25);

# ╔═╡ f484f26e-5cfc-4b6d-9978-f0a844ddedb9
begin
	_background_ring = create_cylinder(mpr, centers_a, centers_b, 12, -25);
	background_ring = Bool.(_background_ring .- cylinder)
end;

# ╔═╡ eae6cdb3-9d6a-4c7a-8d64-0f7a0832ae2c
@bind z PlutoUI.Slider(axes(mpr, 3), default=div(size(mpr, 3), 2), show_value=true)

# ╔═╡ b5b42c05-cc36-4c17-b64f-1fed5a83a6ef
let
	idxs = getindex.(findall(isone, cylinder[:, :, z]), [1 2])
	idxs_ring = getindex.(findall(isone, background_ring[:, :, z]), [1 2])
	α = 0.25

	f = Figure()

	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(transpose(mpr[:, :, z]); colormap = :grays)
	scatter!(idxs[:, 2], idxs[:, 1]; markersize = 2, color = (:red, α), label = "inserts")
	scatter!(idxs_ring[:, 2], idxs_ring[:, 1]; markersize = 3, color = (:blue, α), label = "background")
	axislegend(ax)

	f
end

# ╔═╡ 9b66d146-611e-4433-8426-b8f39815d138
md"""
## Segment Calibration Insert
"""

# ╔═╡ ad13b5a6-1232-48ed-b91e-8b23896a0bd3
begin
	binary_calibration = falses(size(mpr))
	binary_calibration[centers_a...] = true
	binary_calibration = dilate(binary_calibration)
end;

# ╔═╡ 927a753b-e8db-465a-883e-c939e767bede
let
	idxs = getindex.(findall(isone, Int.(binary_calibration[:, :, 160])), [1 2])	
	f = Figure()

	ax = CairoMakie.Axis(f[1, 1])
	heatmap!(transpose(mpr[:, :, 160]); colormap = :grays)
	scatter!(idxs[:, 2], idxs[:, 1]; markersize = 1, color = :red)

	f
end

# ╔═╡ 017cfb8e-6586-4f93-b817-6afab3b10d30
# function _percentage_calcium(voxel_intensity, hu_calcium, hu_heart_tissue)
#     return precentage_calcium = (voxel_intensity - hu_heart_tissue) / (hu_calcium - hu_heart_tissue)
# end

# ╔═╡ 66d46419-3315-4e2f-b11c-71f38ccd1087
# function score(vol::AbstractArray, hu_calcium, hu_heart_tissue, voxel_size, density_calcium)
#     number_calcium_voxels = []
#     for i in axes(vol, 1)
#         for j in axes(vol, 2)
#             for k in axes(vol, 3)
#                 m = vol[i, j, k]
#                 percent = _percentage_calcium(m, hu_calcium, hu_heart_tissue)
#                 push!(number_calcium_voxels, percent)
#             end
#         end
#     end
#     return sum(number_calcium_voxels) * voxel_size * density_calcium
# end

# ╔═╡ d5b62cb5-15ec-421c-a9ba-379b801cde72
md"""
## Remove Outliers (Air)
"""

# ╔═╡ 479ea014-1200-487f-a851-ee3e45339620
function remove_outliers(vector)
    Q1 = quantile(vector, 0.25)
    Q3 = quantile(vector, 0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    return [x for x in vector if x > lower_bound]
end

# ╔═╡ 1a196009-3f5f-45d7-9d3f-ad9da9afcd1b
mpr_clean = remove_outliers(mpr[cylinder]);

# ╔═╡ 71f39103-4ae2-4673-9fe9-7a1a02002d7f
let
	f = Figure(resolution = (1000, 1200))
	ax = CairoMakie.Axis(f[1, 1], title = "Original")
	hist!(mpr[cylinder])

	ax = CairoMakie.Axis(f[2, 1], title = "Clean")
	hist!(mpr_clean)

	f
end

# ╔═╡ 7caa6404-7fdb-4b6e-a976-86f63cab282a
md"""
# Score
"""

# ╔═╡ e28e7c8a-51fd-43d6-9291-05501922db07
md"""
## Ground Truth
"""

# ╔═╡ a8835581-db72-43d3-913d-b62a23607cdc
begin
	gt_density = 0.050 # mg/mm^3

	# π * (diameter/2)^2 * length
	gt_volume = π * (5/2)^2 * 7 # mm3
	gt_mass = gt_density * gt_volume
end

# ╔═╡ a7c997e6-2dcf-487e-8b9d-3524c7a0ea27
md"""
## Volume Fraction
"""

# ╔═╡ 611528e5-d64c-4d01-ab42-e32a8dcda37a
hu_calcium_400 = mean(mpr[binary_calibration])

# ╔═╡ 976ddafe-2806-4ba8-b6ed-32e36939f2c5
ρ_calcium_400 = 0.400 # mg/mm^3

# ╔═╡ 96f14096-400a-43eb-8c24-a8aa393d105a
voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

# ╔═╡ faa8104f-b316-45fa-a733-82611099ea71
hu_heart_tissue_bkg = mean(mpr[background_ring])

# ╔═╡ 729d6ef4-5427-4a1c-af66-460e83505187
vf_mass = score(mpr_clean, hu_calcium_400, hu_heart_tissue_bkg, voxel_size, ρ_calcium_400, VolumeFraction())

# ╔═╡ 88c88bc1-51f3-4c5c-8d5d-5eb9c0e6b64e
md"""
## Agatston
"""

# ╔═╡ f446cc59-1fa9-4747-84bb-59597c4d6039
overlayed_mask = create_mask(mpr, cylinder);

# ╔═╡ ba9b8451-023c-400c-8555-11324b092cbb
kV = header[tag"KVP"]

# ╔═╡ 8f32acdc-3ade-44b9-82f4-bf5031c0a9da
mass_cal_factor = ρ_calcium_400 / hu_calcium_400

# ╔═╡ 0af83bc7-9bfd-4efd-8a6c-3b85eae45e75
a_agatston, a_volume, a_mass = score(overlayed_mask, pixel_size, mass_cal_factor, Agatston(); kV=kV)

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
# ╠═4edf70c6-a181-494d-b996-8bb8a44f1f32
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
# ╠═4b50a86e-fb36-4d26-9460-111f0465dce7
# ╠═6edfa2b7-63ae-4315-ac9f-caf96f0f51ed
# ╠═32d8eec0-892e-455f-9815-26eccc4dd3d6
# ╠═d3282a5e-1a7a-4537-b7f1-cc4a7c99c464
# ╠═24c102fd-1fc0-4c31-9ad1-d9af691159b0
# ╟─8f6a2277-504d-4009-ac15-cca75b1e6575
# ╠═b53a9ab0-38e3-4144-990a-55f1f5cdd33d
# ╠═aa8d9f5a-c272-4dba-894c-653b92501bb3
# ╟─200c67b9-bff2-456c-9b74-0ca1c31ee964
# ╠═97739511-af5a-46f9-a251-09c1e13ce758
# ╠═216c4a93-4acc-4c3e-9856-40ed1daf4611
# ╟─6483fd22-acdb-4325-8d67-bd54c08c9a27
# ╠═649f943d-8838-4258-a108-4de678561d83
# ╠═2e72c9bb-e41c-493f-b62f-9cc60531382a
# ╠═6b1ee364-7611-4fde-9354-aeee659a127f
# ╠═f484f26e-5cfc-4b6d-9978-f0a844ddedb9
# ╟─eae6cdb3-9d6a-4c7a-8d64-0f7a0832ae2c
# ╟─b5b42c05-cc36-4c17-b64f-1fed5a83a6ef
# ╟─9b66d146-611e-4433-8426-b8f39815d138
# ╠═ad13b5a6-1232-48ed-b91e-8b23896a0bd3
# ╟─927a753b-e8db-465a-883e-c939e767bede
# ╠═017cfb8e-6586-4f93-b817-6afab3b10d30
# ╠═66d46419-3315-4e2f-b11c-71f38ccd1087
# ╟─d5b62cb5-15ec-421c-a9ba-379b801cde72
# ╠═479ea014-1200-487f-a851-ee3e45339620
# ╠═1a196009-3f5f-45d7-9d3f-ad9da9afcd1b
# ╟─71f39103-4ae2-4673-9fe9-7a1a02002d7f
# ╟─7caa6404-7fdb-4b6e-a976-86f63cab282a
# ╟─e28e7c8a-51fd-43d6-9291-05501922db07
# ╠═a8835581-db72-43d3-913d-b62a23607cdc
# ╟─a7c997e6-2dcf-487e-8b9d-3524c7a0ea27
# ╠═b6f3279a-1c31-4c42-bac2-1adc9ec317b4
# ╠═611528e5-d64c-4d01-ab42-e32a8dcda37a
# ╠═976ddafe-2806-4ba8-b6ed-32e36939f2c5
# ╠═96f14096-400a-43eb-8c24-a8aa393d105a
# ╠═faa8104f-b316-45fa-a733-82611099ea71
# ╠═729d6ef4-5427-4a1c-af66-460e83505187
# ╟─88c88bc1-51f3-4c5c-8d5d-5eb9c0e6b64e
# ╠═2d0f900a-85e2-4956-b08e-9f85e1d68a60
# ╠═f446cc59-1fa9-4747-84bb-59597c4d6039
# ╠═ba9b8451-023c-400c-8555-11324b092cbb
# ╠═8f32acdc-3ade-44b9-82f4-bf5031c0a9da
# ╠═0af83bc7-9bfd-4efd-8a6c-3b85eae45e75
