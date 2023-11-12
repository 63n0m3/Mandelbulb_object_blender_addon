# Mandelbulb_object_blender_addon
Addon for blender for creation of 3d fractals.
In main is a current working version. All functions should work, however OpenCL crashes my GPU driver on bigger Mandelbulb resolutions.
There are instructions in files how to run it if you have OpenCL developing enviroment or it is possible to just run it in c++ single threaded or just python. I do not provide DLL as of now it has 2 issues.
Version 1.0 is a fully Blender python addon version. Same Mandelbubs are available (all from wikipedia) in 1.0, less mesh cleaning algorithms. Python version is 5-15 times slower than single threaded c++ dll version, Both are incomparable to OpenCL on heavier workloads.
