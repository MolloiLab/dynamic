"""
    offset_interpolate(itp, plane, i, j, k)

Interpolate the intensity values at a new position offset from the original plane.

# Arguments
- `itp::Interpolation`: Interpolation object used for interpolating the intensity values.
- `plane::Meshes.Plane`: Plane object defining the plane in which the MPR is created.
- `i::Int`, `j::Int`, `k::Int`: Indices defining the position of the point in the 3D space.

# Returns
- `Float`: Interpolated intensity value at the new position.
"""
function offset_interpolate(itp, plane, i, j, k)
    itp((plane(i, j) + Meshes.normal(plane)*k).coords...)
end

"""
    create_mpr(centroids, points, r1, r2)

Create a Multi-Planar Reformat (MPR) of a 3D image data.

# Arguments
- `centroids::Array`: Array of centroids defining the center of the MPR.
- `points::Array`: Array of points defining the plane in which the MPR is created.
- `r1::Int`: Radius of the MPR in the i and j dimensions.
- `r2::Int`: Radius of the MPR in the k dimension.

# Returns
- `Array`: The created MPR as a 3D array.
"""
function create_mpr(centroids, points, r1, r2)
    offset = Meshes.Point(centroids) - Meshes.Point(first(points))
    plane = Meshes.Plane((Meshes.Point.(points) .+ Ref(offset))...)
    mpr_itp = linear_interpolation(axes(dcm_heart), dcm_heart; extrapolation_bc = 0.0);
    return [offset_interpolate(mpr_itp, plane, i, j, k) for i in -r1:r1, j in -r1:r1, k in -r2:r2]
end


mpr = create_mpr(centroids, pts, size(dcm_arr, 2)/2, size(dcm_arr, 3)/2);