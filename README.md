# Software-Rasterizer

This is my custom implementation of a software rasterizer, which renders pixels to a pixel buffer and then displays the result to the screen using SDL3.
The project includes my own implementations of vector, matrix, and quaternion math, all provided through custom, lightweight header files.

This project began as a (C++) clone of [Sebastian Lague’s](https://www.youtube.com/@SebastianLague) software rasterizer, but over time it has evolved with enough additions and changes to become its own standalone project.

## Pictures / Showcase

Below are some screenshots and renders captured directly from the software rasterizer,
showcasing the lighting, texture, and alpha effects.

<p align="center">
  <img src="https://github.com/user-attachments/assets/f6d76aa7-dfae-4aec-9776-83cbc7a5cdd7" width="320" alt="LitShaderMonkey">
  <img src="https://github.com/user-attachments/assets/f09819dc-6995-4f64-a82d-f559e7c6f306" width="320" alt="CapriCube">
  <img src="https://github.com/user-attachments/assets/b54f9482-669f-4f74-a8c4-13b0517e8e3e" width="320" alt="CubeWithTransparency">
</p>

## Technical Information

The rasterizer supports loading OBJ and MTL files, which are converted into an internal Object structure.

### Objects

#### Each Object contains:

- Faces
- Transform data: position, rotation, and scale (Vector3)
- Material list
- Shader reference

#### The Face structure holds:

- Local Vertex positions
- Texture coordinates
- Specular map references
- Normal vectors

#### The Material structure holds:

- Base color (Diffuse color)
- Texture data (if applicable)
- Normals
- Opacity
- Other standard OBJ material properties

### Shaders

There are five types of shaders that can be applied to an object:

#### DiffuseShader

A simple shader that returns only the object’s color.
- No lighting
- No texture

#### LitDiffuseShader

Similar to DiffuseShader, but takes into account light objects and calculates color based on the diffuse color and the light’s color and distance.
- Colored light
- Distance based illuminosity
- Supports shadow checks by testing vertex–light intersections (expensive — can significantly reduce FPS)

#### DiffuseAlphaShader

Like DiffuseShader, but supports alpha transparency.
- Blends the object’s color with the background based on opacity
- Rendered after opaque objects to ensure correct blending and depth behavior

#### TextureShader

Uses a texture instead of a diffuse color to color the object’s faces.

#### LitTextureShader

Like TextureShader, but with lighting support (see LitDiffuseShader for details).

### Camera

The camera stores:
- Position
- Rotation (Quaternion) to avoid gimbal lock
- Field of view (FOV)

### Lights

There are three types of lights

#### PointLight

Emits light in all directions from a single point.
- Has a position, color, and intensity

#### DirectionalLight

Emits light in a specific direction, simulating a very distant light source (e.g., sunlight).
- Has a direction and color
- Light rays are considered parallel (infinite range)
- Can impact performance, as generally more vertices need to be checked for shadows compared to a point light

#### Spotlight

Emits light in a cone from a point in space.
- Has a position, direction, color, cone angle, and intensity
- Useful for effects like flashlights or stage lighting

## Render pipeline

The rasterizer processes and displays each frame through the following steps:

### 1. Per-Object Processing

For each object in the scene:

#### 1.1 Matrix Preparation

- Compute the model matrix from position, rotation, and scale
- Get the view matrix from the active camera
- Derive a normal matrix for transforming normals

#### 1.2 Face Processing

For each face in the object:
- Triangulate polygons into triangles using the Ear Clipping algorithm
- Transform vertices from local space -> world space -> view space using the matrices calculated in step 1.1
- Transform normals to world space and normalize them

#### 1.3 Triangle Clipping

- Slice triangles partially behind the near plane of the camera
- Separate opaque and transparent triangles for correct render order (opaque renders first)

### 2. Rasterization

For each triangle from the previous step 1:

#### 2.1 Screen Projection

- Project vertices from view space to screen coordinates, adding size variation based on distance and depth distortion from the FOV.
- (In practice, this projection occurs inside triangleSlicer, but it is logically a separate stage.)

#### 2.2 Backface and Degenerate Checks

- Skip triangles with zero area
- Skip triangles entirely behind the camera

#### 2.3 Bounding Box & Edge Rules

- Compute the screen-space bounding box of the triangle
- Apply the top-left edge bias rule to avoid pixel overlap between adjacent triangles

#### 2.4 Pixel Iteration

- Loop over each pixel inside the bounding box
- Use barycentric coordinates to determine pixel coverage and interpolate attributes

### 3. Depth Testing

- Compare interpolated pixel depth against the depth buffer
- Discard the pixel if something closer has already been drawn (reason opaque geometry is rendered first).

### 4. Shading

- Interpolate texture coordinates, normals, and world-space position per pixel.
- Call the object’s shader to determine the final pixel color.
- The shader receives:
  - The material
  - Interpolated attributes
  - Light list
  - All scene triangles (for shadow checks in lit shaders)
- Example:
  - DiffuseShader → returns the material’s diffuse color.
  - Lit shaders → perform per-pixel lighting and shadow checks before blending with diffuse or texture color.

### 5. Writing to Buffers

If the shader outputs a visible color:
- Write it to the draw buffer
- Update the depth buffer with the new pixel depth
