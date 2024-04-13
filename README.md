# dynamic

## Data Structure
### Folders
- `a`: 1.2 mm diameter, 50 mg/cc density
- `b`: 3.0 mm diameter, 50 mg/cc density
- `c-4x`: 5.0 mm diameter, 50 mg/cc density
- `d`: 1.2 mm diameter, 100 mg/cc density
- `e`: 3.0 mm diameter, 100 mg/cc density
- `f`: 5.0 mm diameter, 100 mg/cc density
- `h`: 1.2 mm diameter, 250 mg/cc density
- `i`: 3.0 mm diameter, 250 mg/cc density
- `k`: 5.0 mm diameter, 250 mg/cc density
- `l`: 1.2 mm diameter, 400 mg/cc density
- `m`: 3.0 mm diameter, 400 mg/cc density
- `n`: 5.0 mm diameter, 400 mg/cc density

### Phantom Insert Layout
![alt text](images/phantom_layout.png "Phantom Insert Layout")

Notice how the phantom is set up in a way for easy automatic segmentation. Specifically the phantom contains mostly background inserts and only one row of inserts for measurements, without outer inserts of the highest density large inserts (`n`) that are not meant for measurement but specifically just for detection of automatic segmentation. Then it contains another set up background inserts, then three measurement inserts with varying amounts of calcium (`a` - `n`)

*NOTE: In order to avoid the effect of air as much as possible, the inserts were surrounded by silly putty that is specifically around the hounsfield unit (HU)/intensity of the CIRS background tissue (~ 50 HU). This is likely imperfect and will have some effect on the results, but it is much better that having air (~ -1000 HU)*

### More information
The root folders (e.g. a, b, c-4x, d, e, ...) containing the data are all named according to the inserts they contain (see above). Within these folders there are many more subfolders, which I have not yet thoroughly investigated, but correspond to different scan parameters. For example, dual energy scans are indicated via "Spectral ..." and the corresponding monoenergetic reconstruction is indicated via "50 Mono", "60 Mono", etc. Slice thickness reconstructions are indicated via "0.5", "1.0", etc. Somewhere within the folder, cardiac phase reconstruction is indicated via the percentage. We also performed different types of reconstruction algorithms, e.g. filtered back projection, iterative, etc., but these are not indicated within the folders themselves, so this needs to be determined via the DICOM header or some other manner. Potentially, all of this information can be determined via the DICOM header.

The phantom is a dynamic phantom so the scans are all either 0, 60, or 90 beats per minute (bpm) [I THINK]. This is hopefully determined via the header too but I am not positive right now

### General Table of Inserts
#### Not specific to this study but the overall CIRS phantom
| Material                | Densities:  | Diameters:                     |                            |                     |                         |                     |
|-------------------------|-------------|--------------------------------|----------------------------|---------------------|-------------------------|---------------------|
|                         | 50 mg/cc    | 1.2mm Dark pink with a cream   | 3.0mmDark pink with a cream| 5.0mm               |                         |                     |
|                         |             | Ø1.2mm core (code A)           | Ø3 mm core (code B)        | Cream (code G)      |                         |                     |
| Hydroxyapatite          |             | Batch:                         | Qty:                       | Batch:              | Qty: 6                  |                     |
|                         | 100 mg/cc   | 1.2mm Dark pink with a cream   | 3.0mm Dark pink with a cream| 5.0mm Cream (lighter|                         |                     |
| (calcification targets) |             | Ø1.2mm core (code D)           | Ø3 mm core (code E)         | than for 50mg/cc)   |                         |                     |
|                         |             | Batch:                         | Qty:                       | (code F)            |                         |                     |
| Background Batch:       |             |                                |                            | Batch:              | Qty: 4                  |                     |
|                         | 250 mg/cc   | 1.2mm Dark pink with a cream   | 3.0mm Dark pink with a cream| 5.0mm Cream (lighter|                         |                     |
|                         |             | Ø1.2mm core (code H)           | Ø3 mm core (code I)         | than for 100mg/cc)  |                         |                     |
|                         |             | Batch:                         | Qty:                       | (code K)            |                         |                     |
|                         |             |                                |                            | Batch:              | Qty: 4                  |                     |
|                         | 400 mg/cc   | 1.2mm Dark pink with a cream   | 3.0mm Dark pink with a cream| 5.0mm Cream (lighter|                         |                     |
|                         |             | Ø1.2mm core (code L)           | Ø3 mm core (code M)         | than for 250mg/cc)  |                         |                     |
|                         |             | Batch:                         | Qty: 6                     | (code N)            |                         |                     |
|                         |             |                                |                            | Batch:              | Qty: 4                  |                     |
|                         | 0. 5 mg/cc  | 1.2mm Dark pink with a cream   | 3.0mm Dark pink with a cream| 5.0mm Translucent   |                         |                     |
|                         |             | Ø1.2mm core (code O)           | Ø3mm core (code P)          | red                 |                         |                     |
|                         |             | Batch:                         | Qty: 6                     | (code R)            |                         |                     |
| Iodine                  | 1.0 mg/cc   | 1.2mm Dark pink with a cream   | 3.0mm Dark pink with a cream| Batch:              | Qty: 6                  |                     |
|                         |             | Ø1.2mm core (code S)           | Ø3mm core (code T)          | 5.0mm Translucent   |                         |                     |
| Background Batch:       |             | Batch:                         | Qty: 6                     | blue with a green   |                         |                     |
|                         |             |                                |                            | nuance (code U)     |                         |                     |
|                         | 5.0 mg/cc   | 1.2mm Dark pink with a cream   | 3.0mm Dark pink with a cream| Batch:              | Qty: 6                  |                     |
|                         |             | Ø1.2mm core (code V)           | Ø3mm core (code X)          | 5.0mm Translucent   |                         |                     |
|                         |             | Batch:                         | Qty: 3                     | Yellow (code Y)     |                         |                     |
|                         |             |                                |                            | Batch:              | Qty: 6                  |                     |
|                         | 10 mg/cc    | 1.2mm Dark pink with a cream   | 3.0mm Dark pink with a cream| 5.0mm Translucent   |                         |                     |
|                         |             | Ø1.2mm core (code Z)           | Ø3mm core (code 2)          | purple              |                         |                     |
|                         |             | Batch:                         | Qty: 6                     | (code 8)            |                         |                     |
|                         |             |                                |                            | Batch:              | Qty: 6                  |                     |
| Solid Background Insert | Batch       |                                | Qty: 12                    |                     |                         |                     |

## DrWatson.jl stuff
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> dynamic

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "dynamic"
```
which auto-activate the project and enable local path handling from DrWatson.