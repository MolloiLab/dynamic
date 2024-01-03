"""
    load_dcm_array(dcm_data)
Given some DICOM.DICOMData, `load_dcm_array` loads the pixels
of each slice into a 3D array and returns the array
"""
function load_dcm_array(dcm_data)
    return array = cat(
        [dcm_data[i][(0x7fe0, 0x0010)] for i in eachindex(dcm_data)]...; dims=3
    )
end