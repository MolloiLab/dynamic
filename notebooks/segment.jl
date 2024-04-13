### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ c2170328-dc64-4617-8cb7-2cef3d463f02
using Pkg; Pkg.instantiate()

# ╔═╡ 011b0ae2-48d4-4141-8bdb-c94d87ef0a38
using DrWatson

# ╔═╡ c66fda61-be14-4a08-810a-0d75ba617786
# ╠═╡ show_logs = false
@quickactivate "dynamic"

# ╔═╡ 7135f9c6-ba3c-492b-8d1f-f5928003d5ca
using DataFrames: DataFrame, rename!, nrow

# ╔═╡ 49232c37-a857-4946-a0df-af2ce3135229
using StaticArrays: SVector

# ╔═╡ e0fb472b-8006-4b03-85ba-7fda0080f740
using IterTools: product

# ╔═╡ 3f62a1bb-b764-4dd3-a491-86d7bb33acc7
using LinearAlgebra: normalize, dot

# ╔═╡ daa6f5c9-a613-4032-8d75-01883d3c6cc4
using CairoMakie: Figure, Axis, heatmap!, heatmap, scatter!, axislegend, hist!

# ╔═╡ a3b12838-bfb6-4b20-8fcf-b77ea2ae3897
using DICOM: dcmdir_parse

# ╔═╡ 0781d60f-3667-4def-9b7f-54a4ef5d8980
using PlutoUI: TableOfContents, Slider

# ╔═╡ 5169852a-4e29-4e3d-9aa0-4626138c0adf
using StatsBase: countmap, quantile

# ╔═╡ d3b0604e-bb41-4d83-9518-3d51809bde48
using ImageMorphology: component_centroids, label_components, component_boxes, dilate

# ╔═╡ 53884809-3f66-4639-94f7-ed84a3a60883
using Statistics: mean, norm, std

# ╔═╡ fb0d54d2-13b4-46d9-be9a-a0e0285cf5a6
using CalciumScoring: score, VolumeFraction, Agatston

# ╔═╡ 3cb44528-6fcd-419b-8862-3487fced23b9
include(srcdir("dicomutils.jl"));

# ╔═╡ 3ac15a2c-382c-4a24-99fa-9192ab0b5f30
include(srcdir("active_contours.jl"));

# ╔═╡ 2d0f900a-85e2-4956-b08e-9f85e1d68a60
include(srcdir("masks.jl"));

# ╔═╡ 0b0d1327-83a3-4ed0-a5d9-ecbdada1bd64
md"""
# Set Up
"""

# ╔═╡ 446aff6b-3b74-4bb9-b84e-27dc898b90ec
md"""
## Environment & Imports
"""

# ╔═╡ 7e701979-31f6-45f2-a4d6-6928d93ad042
TableOfContents()

# ╔═╡ 984cb619-c918-4558-a0ee-6b82cf03bf76
md"""
## Docs
"""

# ╔═╡ f40642fe-2fc4-4bde-bf9f-690c0c9e30eb
begin
	insert_radii = [1.2, 3.0, 5.0] ./ 2 # mm
	insert_densities = [0.05, 0.1, 0.25, 0.4] # mg/mm^3
	heart_rates = [0, 60, 90] # bpm
	
	scan_types = ["de_80kv", "de_135kv", "se_80kv", "se_120kv", "se_135kv"]
	slice_thicknesses = [0.5, 1.0, 3.0] # mm
	reconstruction_types = ["fbp", "aidr"]
end;

# ╔═╡ 87721c59-de02-45e8-8077-1709ab413696
prod = product(insert_radii, insert_densities, heart_rates, scan_types, slice_thicknesses, reconstruction_types);

# ╔═╡ 95cb337c-8b57-424b-9399-c966dab6b53c
begin
	df = DataFrame(collect(prod))
	rename!(df, [:insert_radii, :insert_densities, :heart_rates, :scan_types, :slice_thickness, :reconstruction_types])
	sort!(df, [:insert_radii, :insert_densities, :heart_rates])
end

# ╔═╡ 9f6c2532-07d1-4b6d-8606-e6557bc26439
total_accessions = nrow(df) / 30

# ╔═╡ 340a1e3d-9d29-4afe-88d3-2e38da413595
md"""
# Load DICOMs
"""

# ╔═╡ 13c0d236-50cf-4c97-9604-13699cdd7979
root = "/Volumes/USB DISK/all/b"

# ╔═╡ b7d678bf-22b1-44db-bdc1-ea13e71c82e0
output_path = joinpath(root, "Cardiac 0.5", "14")

# ╔═╡ b8a51c87-90cf-4259-9061-4e8dc8e0632c
dcms = dcmdir_parse(output_path)

# ╔═╡ b38f92ac-ad86-404b-bf33-118d7f886b36
dcm_arr = load_dcm_array(dcms);

# ╔═╡ 33eae10c-3777-4c53-8245-72932a28069d
header = dcms[1].meta;

# ╔═╡ 231939ec-ee93-4363-a962-ac7517b0bb0a
header

# ╔═╡ 3ad9eb7c-9241-4a06-8308-b6ba33aa492f
md"""
!!! warning
	The `pixel_size` should automatically be extracted from the DICOM header information using `get_pixel_size(header)` (from the src file `dicomutils.jl`). But the dicom tag associated with the pixel size is not in the header for some reason, so right now this is being hardcoded below. This should be investigated further
"""

# ╔═╡ c7cd528e-8ee0-47ee-aaf6-47158532e2ae
pixel_size = [0.5, 0.5, 0.5]

# ╔═╡ 6321ab34-df67-410a-8776-850b4b456058
md"""
## Visualize
"""

# ╔═╡ 3fd05c3f-eb99-4d26-8b55-eb96b1513bc1
@bind a Slider(axes(dcm_arr, 3), default=160, show_value=true)

# ╔═╡ 496773b5-a3bd-4ed7-aba5-7effe8cd958b
let
	f = Figure()

	Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, a], colormap=:grays)

	f
end

# ╔═╡ 01d783da-dfcf-4031-bd21-2dc5a4e014d2
md"""
# Mask Heart Function
"""

# ╔═╡ 63acd515-78aa-4f3f-a472-4c6dfd2fea9d
function erode_mask(img, num_erosions)
	new_img = img
	i = 0
	while i < num_erosions
		new_img = erode(new_img)
		i += 1
	end
	return new_img
end

# ╔═╡ 61201eca-e986-41ad-9d90-2ddfbc4a5a90
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

# ╔═╡ 7effead4-1993-47d2-9bb6-4c58d0d2cc21
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

# ╔═╡ 39a8257c-2664-4779-ba48-6835fbb7b290
begin
	half_x, half_y = size(dcm_arr, 1) ÷ 2, size(dcm_arr, 2) ÷ 2
	init_circle = create_circle_mask(dcm_arr[:, :, 3], (half_x, half_y), 140)

	init_mask = BitArray(undef, size(dcm_arr))
	for z in axes(dcm_arr, 3)
		init_mask[:, :, z] = init_circle
	end

	init_mask = init_mask .* initial_level_set(size(init_mask))
end;

# ╔═╡ 23a3c740-6477-42b4-aa4e-990130b2bf64
heart_cv = chan_vese(dcm_arr; init_level_set = init_mask);

# ╔═╡ ddcf30bb-0db0-4c1c-bec0-c1eec9289bba
centroids = centroids_from_mask(heart_cv)

# ╔═╡ 020c0257-1d58-4f26-aefb-b019024b4629
md"""
!!! warning
	This radius for the heart mask is currently hardcoded to `100` as seen below. This likely wont be good enough and needs to be somewhat more dynamic to work for every type of scan. This should be investigated more
"""

# ╔═╡ d77655b5-2249-4fe7-91c5-5f2bc69f2aba
heart_rad = 100

# ╔═╡ 70461a00-bad3-432d-a2b4-3fe0d5dc0817
heart_mask = create_circle_mask(dcm_arr[:, :, 3], centroids, heart_rad);

# ╔═╡ 33eed811-b9a5-4589-8f6c-721d48ef8c9e
md"""
## Visualize
"""

# ╔═╡ 7773815d-c91d-49a4-a0dd-46128ee6b28b
@bind z2 Slider(axes(dcm_arr, 3), default=130, show_value=true)

# ╔═╡ bbadeeb3-ac1e-4a5a-a70d-69e0c0d57a02
let
	f = Figure()

	Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, z2]; colormap=:grays)
	heatmap!(heart_mask; colormap=(:jet, 0.3))

	f
end

# ╔═╡ 40dee94c-a483-4197-9573-2f69544648ba
md"""
# Inserts
"""

# ╔═╡ 8fba68fa-ff44-4962-9d6f-9b274ba102b1
dcm_heart = dcm_arr .* heart_mask;

# ╔═╡ 6fba9862-e8e5-4faf-8443-ff85c06c30d1
@bind d Slider(axes(dcm_heart, 3), default=130, show_value=true)

# ╔═╡ 97e205d4-9ad5-4cf4-9b64-fd0d09d5bc32
let
	f = Figure()

	Axis(f[1, 1])
	heatmap!(dcm_heart[:, :, d], colormap=:grays)
	
	f
end

# ╔═╡ 200c67b9-bff2-456c-9b74-0ca1c31ee964
md"""
## Segment Calcium Inserts
"""

# ╔═╡ 97739511-af5a-46f9-a251-09c1e13ce758
function get_insert_centers(dcm, threshold)
	dcm_slice_thresh = dcm .> threshold
	
	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(dcm_slice_thresh)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# Find the value with the most occurrences
	most_common_value_a, most_common_value_b = sort(collect(pairs(counts)), by=x->x[2], rev=true)

	# Find the indices of the most common value in the original array
	most_common_indices_a = findall(cc_labels .== most_common_value_a[1])

	# Create boolean array from new cartesian indices
	bool_arr_a = zeros(size(dcm_slice_thresh))
	for i in most_common_indices_a
		bool_arr_a[i] = 1
	end
	centroids_a = Int.(round.(component_centroids(label_components(bool_arr_a))[end]))
	box_a = component_boxes(label_components(bool_arr_a))

	# Find the indices of the most common value in the original array
	most_common_indices_b = findall(cc_labels .== most_common_value_b[1])

	# Create boolean array from new cartesian indices
	bool_arr_b = zeros(size(dcm_slice_thresh))
	for i in most_common_indices_b
		bool_arr_b[i] = 1
	end
	centroids_b = Int.(round.(component_centroids(label_components(bool_arr_b))[end]))

	# centers_a, centers_b = [centroids_a..., z], [centroids_b..., z]
	centers_a, centers_b = centroids_a, centroids_b
	return centers_a, centers_b
	
end

# ╔═╡ 216c4a93-4acc-4c3e-9856-40ed1daf4611
centers_a, centers_b = get_insert_centers(dcm_heart, 200);

# ╔═╡ 0d5cfc80-8b9a-48bb-ba39-a44a0be16a10
centers_a

# ╔═╡ 15775655-9f2a-4d56-bf66-e919a7dc4284
centers_b

# ╔═╡ 852bf7a3-04f6-47a0-9009-dff9f186f9cf
@bind z3 Slider([centers_a[3], centers_b[3]], show_value = true)

# ╔═╡ 6483fd22-acdb-4325-8d67-bd54c08c9a27
let
	msize = 10
	f = Figure()

	ax = Axis(f[1, 1])
	# z = div(size(dcm_heart, 3), 2)
	# heatmap!(mpr[:, :, z]; colormap=:grays)
	heatmap!(dcm_heart[:, :, z3]; colormap=:grays)
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
cylinder = create_cylinder(dcm_heart, centers_a, centers_b, 8, -25);

# ╔═╡ f484f26e-5cfc-4b6d-9978-f0a844ddedb9
begin
	_background_ring = create_cylinder(dcm_heart, centers_a, centers_b, 12, -25);
	background_ring = Bool.(_background_ring .- cylinder)
end;

# ╔═╡ eae6cdb3-9d6a-4c7a-8d64-0f7a0832ae2c
@bind z Slider(axes(dcm_heart, 3), default=div(size(dcm_heart, 3), 2), show_value=true)

# ╔═╡ b5b42c05-cc36-4c17-b64f-1fed5a83a6ef
let
	idxs = getindex.(findall(isone, cylinder[:, :, z]), [1 2])
	idxs_ring = getindex.(findall(isone, background_ring[:, :, z]), [1 2])
	α = 0.25

	f = Figure()

	ax = Axis(f[1, 1])
	heatmap!(transpose(dcm_arr[:, :, z]); colormap = :grays)
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
	binary_calibration = falses(size(dcm_heart))
	binary_calibration[centers_a...] = true
	binary_calibration = dilate(binary_calibration)
end;

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
dcm_heart_clean = remove_outliers(dcm_heart[cylinder]);

# ╔═╡ 71f39103-4ae2-4673-9fe9-7a1a02002d7f
let
	f = Figure(size = (1000, 1200))
	ax = Axis(f[1, 1], title = "Original")
	hist!(dcm_heart[cylinder])

	ax = Axis(f[2, 1], title = "Clean")
	hist!(dcm_heart_clean)

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

# ╔═╡ 42b57843-fd04-4671-9257-c17899d59021
md"""
!!! warning
	The ground truth data is dependent on the scan that is being uploaded. This means that the ground truth needs to be dynamic, based on the folder (e.g. phantom scans A, or B, or ...), but right now it is hard coded and likely not accurate depending on the folder you are in

	*TODO: Update this to be dynamic/folder specific*
"""

# ╔═╡ a8835581-db72-43d3-913d-b62a23607cdc
begin
	gt_density = 0.10 # mg/mm^3

	# π * (diameter/2)^2 * length
	gt_volume = π * (5/2)^2 * 7 # mm3
	gt_mass = gt_density * gt_volume
end

# ╔═╡ a7c997e6-2dcf-487e-8b9d-3524c7a0ea27
md"""
## Volume Fraction
"""

# ╔═╡ 611528e5-d64c-4d01-ab42-e32a8dcda37a
hu_calcium_400 = mean(dcm_heart[binary_calibration])

# ╔═╡ a81a5f9d-359b-454d-bc29-d8770bb1a078
std(dcm_heart[binary_calibration])

# ╔═╡ 976ddafe-2806-4ba8-b6ed-32e36939f2c5
ρ_calcium_400 = 0.400 # mg/mm^3

# ╔═╡ 96f14096-400a-43eb-8c24-a8aa393d105a
voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

# ╔═╡ faa8104f-b316-45fa-a733-82611099ea71
hu_heart_tissue_bkg = mean(dcm_heart[background_ring])

# ╔═╡ 729d6ef4-5427-4a1c-af66-460e83505187
vf_mass = score(dcm_heart_clean, hu_calcium_400, hu_heart_tissue_bkg, voxel_size, ρ_calcium_400, VolumeFraction())

# ╔═╡ 88c88bc1-51f3-4c5c-8d5d-5eb9c0e6b64e
md"""
## Agatston
"""

# ╔═╡ 3b2a39ca-b19a-4333-86b5-2672bced55d1
begin
	overlayed_mask = zeros(size(dcm_arr))
	for z in axes(dcm_arr, 3)
		overlayed_mask = dcm_arr[:, :, z] .* heart_mask
	end
end

# ╔═╡ ba9b8451-023c-400c-8555-11324b092cbb
kV = 120

# ╔═╡ 8f32acdc-3ade-44b9-82f4-bf5031c0a9da
mass_cal_factor = ρ_calcium_400 / hu_calcium_400

# ╔═╡ 0af83bc7-9bfd-4efd-8a6c-3b85eae45e75
a_agatston, a_volume, a_mass = score(overlayed_mask, pixel_size, mass_cal_factor, Agatston(); kV=kV)

# ╔═╡ Cell order:
# ╟─0b0d1327-83a3-4ed0-a5d9-ecbdada1bd64
# ╟─446aff6b-3b74-4bb9-b84e-27dc898b90ec
# ╠═011b0ae2-48d4-4141-8bdb-c94d87ef0a38
# ╠═c66fda61-be14-4a08-810a-0d75ba617786
# ╠═c2170328-dc64-4617-8cb7-2cef3d463f02
# ╠═7135f9c6-ba3c-492b-8d1f-f5928003d5ca
# ╠═49232c37-a857-4946-a0df-af2ce3135229
# ╠═e0fb472b-8006-4b03-85ba-7fda0080f740
# ╠═3f62a1bb-b764-4dd3-a491-86d7bb33acc7
# ╠═daa6f5c9-a613-4032-8d75-01883d3c6cc4
# ╠═a3b12838-bfb6-4b20-8fcf-b77ea2ae3897
# ╠═0781d60f-3667-4def-9b7f-54a4ef5d8980
# ╠═5169852a-4e29-4e3d-9aa0-4626138c0adf
# ╠═d3b0604e-bb41-4d83-9518-3d51809bde48
# ╠═53884809-3f66-4639-94f7-ed84a3a60883
# ╠═fb0d54d2-13b4-46d9-be9a-a0e0285cf5a6
# ╠═3cb44528-6fcd-419b-8862-3487fced23b9
# ╠═3ac15a2c-382c-4a24-99fa-9192ab0b5f30
# ╠═2d0f900a-85e2-4956-b08e-9f85e1d68a60
# ╠═7e701979-31f6-45f2-a4d6-6928d93ad042
# ╟─984cb619-c918-4558-a0ee-6b82cf03bf76
# ╠═f40642fe-2fc4-4bde-bf9f-690c0c9e30eb
# ╠═87721c59-de02-45e8-8077-1709ab413696
# ╠═95cb337c-8b57-424b-9399-c966dab6b53c
# ╠═9f6c2532-07d1-4b6d-8606-e6557bc26439
# ╟─340a1e3d-9d29-4afe-88d3-2e38da413595
# ╠═13c0d236-50cf-4c97-9604-13699cdd7979
# ╠═b7d678bf-22b1-44db-bdc1-ea13e71c82e0
# ╠═b8a51c87-90cf-4259-9061-4e8dc8e0632c
# ╠═b38f92ac-ad86-404b-bf33-118d7f886b36
# ╠═33eae10c-3777-4c53-8245-72932a28069d
# ╠═231939ec-ee93-4363-a962-ac7517b0bb0a
# ╟─3ad9eb7c-9241-4a06-8308-b6ba33aa492f
# ╠═c7cd528e-8ee0-47ee-aaf6-47158532e2ae
# ╟─6321ab34-df67-410a-8776-850b4b456058
# ╟─3fd05c3f-eb99-4d26-8b55-eb96b1513bc1
# ╟─496773b5-a3bd-4ed7-aba5-7effe8cd958b
# ╟─01d783da-dfcf-4031-bd21-2dc5a4e014d2
# ╠═63acd515-78aa-4f3f-a472-4c6dfd2fea9d
# ╠═61201eca-e986-41ad-9d90-2ddfbc4a5a90
# ╠═7effead4-1993-47d2-9bb6-4c58d0d2cc21
# ╠═39a8257c-2664-4779-ba48-6835fbb7b290
# ╠═23a3c740-6477-42b4-aa4e-990130b2bf64
# ╠═ddcf30bb-0db0-4c1c-bec0-c1eec9289bba
# ╟─020c0257-1d58-4f26-aefb-b019024b4629
# ╠═d77655b5-2249-4fe7-91c5-5f2bc69f2aba
# ╠═70461a00-bad3-432d-a2b4-3fe0d5dc0817
# ╟─33eed811-b9a5-4589-8f6c-721d48ef8c9e
# ╟─7773815d-c91d-49a4-a0dd-46128ee6b28b
# ╟─bbadeeb3-ac1e-4a5a-a70d-69e0c0d57a02
# ╟─40dee94c-a483-4197-9573-2f69544648ba
# ╠═8fba68fa-ff44-4962-9d6f-9b274ba102b1
# ╟─6fba9862-e8e5-4faf-8443-ff85c06c30d1
# ╟─97e205d4-9ad5-4cf4-9b64-fd0d09d5bc32
# ╟─200c67b9-bff2-456c-9b74-0ca1c31ee964
# ╠═97739511-af5a-46f9-a251-09c1e13ce758
# ╠═216c4a93-4acc-4c3e-9856-40ed1daf4611
# ╠═0d5cfc80-8b9a-48bb-ba39-a44a0be16a10
# ╠═15775655-9f2a-4d56-bf66-e919a7dc4284
# ╟─852bf7a3-04f6-47a0-9009-dff9f186f9cf
# ╟─6483fd22-acdb-4325-8d67-bd54c08c9a27
# ╠═649f943d-8838-4258-a108-4de678561d83
# ╠═2e72c9bb-e41c-493f-b62f-9cc60531382a
# ╠═6b1ee364-7611-4fde-9354-aeee659a127f
# ╠═f484f26e-5cfc-4b6d-9978-f0a844ddedb9
# ╟─eae6cdb3-9d6a-4c7a-8d64-0f7a0832ae2c
# ╟─b5b42c05-cc36-4c17-b64f-1fed5a83a6ef
# ╟─9b66d146-611e-4433-8426-b8f39815d138
# ╠═ad13b5a6-1232-48ed-b91e-8b23896a0bd3
# ╟─d5b62cb5-15ec-421c-a9ba-379b801cde72
# ╠═479ea014-1200-487f-a851-ee3e45339620
# ╠═1a196009-3f5f-45d7-9d3f-ad9da9afcd1b
# ╟─71f39103-4ae2-4673-9fe9-7a1a02002d7f
# ╟─7caa6404-7fdb-4b6e-a976-86f63cab282a
# ╟─e28e7c8a-51fd-43d6-9291-05501922db07
# ╟─42b57843-fd04-4671-9257-c17899d59021
# ╠═a8835581-db72-43d3-913d-b62a23607cdc
# ╟─a7c997e6-2dcf-487e-8b9d-3524c7a0ea27
# ╠═611528e5-d64c-4d01-ab42-e32a8dcda37a
# ╠═a81a5f9d-359b-454d-bc29-d8770bb1a078
# ╠═976ddafe-2806-4ba8-b6ed-32e36939f2c5
# ╠═96f14096-400a-43eb-8c24-a8aa393d105a
# ╠═faa8104f-b316-45fa-a733-82611099ea71
# ╠═729d6ef4-5427-4a1c-af66-460e83505187
# ╟─88c88bc1-51f3-4c5c-8d5d-5eb9c0e6b64e
# ╠═3b2a39ca-b19a-4333-86b5-2672bced55d1
# ╠═ba9b8451-023c-400c-8555-11324b092cbb
# ╠═8f32acdc-3ade-44b9-82f4-bf5031c0a9da
# ╠═0af83bc7-9bfd-4efd-8a6c-3b85eae45e75
