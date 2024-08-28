# 3D Isotropic Delaunay Mesh Generation (Version 0.3)
 
Chen Jianjun  Center for Engineering & Scientific Computation,
Zhejiang University, P. R. China
Copyright reserved, 2005, 10, 26

For further information, please conctact
- Tel: +86-571-87953165
- Fax: +86-571-87953167

Mail: zdchenjj@yahoo.com.cn

## Build on Unix

First clone the source code to your computer, and type
```Bash
mkdir build
cd build
cmake ..
make
```
to build the project. 

## Build on windows
### x64
For example:
```
cmake -G"Visual Studio 11 2012 Win64"
```

### win32
For example:
```
cmake -G"Visual Studio 11 2012"
```

## Run demo
Built the project, you will get an executable named `DTIso3D`,
please copy the mesh files to the build directory, and type
```Bash
./DTIso3D <proj>
```
to execute it, here `<proj>` means the common prefices of the mesh files.
