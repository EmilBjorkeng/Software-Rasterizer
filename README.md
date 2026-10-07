# Software-Rasterizer

This is my custom implementation of a software rasterizer, which renders pixels to a pixel buffer and then displays the result to the screen using SDL3.
The project includes my own implementations of vector, matrix, and quaternion math, all provided through custom, lightweight header files.

This project began as a (C++) clone of [Sebastian Lague’s](https://www.youtube.com/@SebastianLague) software rasterizer, but over time it has evolved with enough additions and changes to become its own standalone project.

## Pictures / Showcase

Below are some screenshots and renders captured directly from the software rasterizer,
showcasing the lighting, texture, and alpha effects.

<p align="center">
  <img src="https://github.com/user-attachments/assets/f6d76aa7-dfae-4aec-9776-83cbc7a5cdd7" width="300" alt="LitShaderMonkey">
  <img src="https://github.com/user-attachments/assets/f09819dc-6995-4f64-a82d-f559e7c6f306" width="300" alt="CapriCube">
  <img src="https://github.com/user-attachments/assets/b54f9482-669f-4f74-a8c4-13b0517e8e3e" width="300" alt="CubeWithTransparency">
</p>

### Installing dependencies

**Arch Linux**
```sh
sudo pacman -S --needed base-devel pkgconf sdl3 sdl3_ttf
```

**Debian / Ubuntu** (Ubuntu 25.10+ or another release that packages SDL3)
```sh
sudo apt install build-essential pkg-config libsdl3-dev libsdl3-ttf-dev
```

**Fedora**
```sh
sudo dnf install gcc-c++ make pkgconf-pkg-config SDL3-devel SDL3_ttf-devel
```

### Building and running

```sh
git clone https://github.com/EmilBjorkeng/Software-Rasterizer.git
cd Software-Rasterizer
make run
```

`make` builds without running, `make debug` builds with debug symbols,
and `make clean` removes build output.

> **Note:** Run the program from the repository root. Assets (models, textures, font)
> are loaded with relative paths like `assets/Cube.obj`.

