# vmesh build workflow

This repository supports two vmesh usage modes:

1. Build `vmesh_test.exe` locally and run API-level test input files.
2. Build the `vmesh` static library in the cloud desktop project and let the outer program call APIs from `vmesh/include/vmesh.h`.

## Required third-party headers

The code includes these header-only dependencies:

- `Eigen`: expected at `third/eigen` locally, or `extern/eigen` in the outer project.
- `libigl`: expected at `third/libigl/include` locally, or `extern/libigl/include` in the outer project.
- `spdlog`: expected at `third/spdlog/include` locally, or `extern/spdlog/include` in the outer project.

The current local checkout already has `third/eigen`, but does not have `third/libigl` or `third/spdlog`.

## Local executable test

Configure:

```powershell
cmake -S . -B build
```

Build only the vmesh test executable:

```powershell
cmake --build build --config Release --target vmesh_main
```

The executable target is `vmesh_main`, and its output file name is `vmesh_test.exe`.

Run it with a debug input emitted by the API:

```powershell
.\build\vmesh\Release\vmesh_test.exe BoundaryMeshing.txt boundary_layer.vtk top_layer.vtk
```

If output paths are omitted, it uses:

- `BoundaryMeshing.txt`
- `boundary_layer.vtk`
- `top_layer.vtk`

## Cloud desktop library build

When this repo is included by the outer project, keep using the `vmesh` target for the API library:

```cmake
add_subdirectory(path/to/blmesh)
target_link_libraries(outer_program PRIVATE vmesh)
target_include_directories(outer_program PRIVATE path/to/blmesh/vmesh/include)
```

If the outer project expects one merged static library instead of many linked targets, use the existing `getLib.bat` idea after building Release/Debug. It merges `vmesh` and its dependent static libraries into `lib/vmesh_r.lib` or `lib/vmesh_d.lib`.

## Optional target switch

The local executable is controlled by:

```powershell
cmake -S . -B build -DVMESH_BUILD_EXE=ON
```

For cloud builds that only need the library:

```powershell
cmake -S . -B build -DVMESH_BUILD_EXE=OFF
```
