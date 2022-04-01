# Example-based stylization of 3D renderings
INF584 project - Simple path tracer, rasterizer, and example-based stylization of 3D renderings

This project features a rasterizer, a global illumination algorithm, as well as an implementation of the example-based synthesis algorithm  _Algorithm 1._ described in section _3.3 Synthesis algorithm_, of the paper [Stylit](https://dcgi.fel.cvut.cz/home/sykorad/stylit.html).
The rendering algorithm is a Montecarlo pathtracer, based on source code of a raytracer provided during the INF584 - Image Synthesis course by Tamy Boubekeur. 

# Instructions

**Compiling**

```bash
cd <path-to-MyRenderer>
mkdir build
cd build
cmake ..
cd ..
cmake --build build -–config Release 
```

The main program and its dependencies will be generated using your local compiler.

**Running**
```bash
cd <path-to-MyRenderer>
./MyRenderer <mesh-file.off>
```
Example Meshes in .off format are located in `Myrenderer/Resources/Models`.


# Libraries (included)
- OpenGL, for accessing your graphics processor
- GLAD, for accessing modern OpenGL extensions
- GLFW to interface OpenGL and the window system of your operating system
- GLM, for the basic mathematical tools (vectors, matrices, etc.).
- Eigen, for matrices

# Example Results
