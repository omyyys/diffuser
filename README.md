# diffuser

**diffuser** is an interactive 2D fluid simulator written in C++ using [raylib](https://github.com/raysan5/raylib) and [raygui](https://github.com/raysan5/raygui). Drag the mouse across the canvas to add fluid and set it in motion, then experiment with viscosity, diffusion and colour. The simulation can display density, velocity or vorticity.

## Building

The project uses **Clang**, **Make** and **C++20**. Raylib 5.5 and raygui are included in the repository.

On Linux, you can use the bundled libraries in `raylib/lib`:

```sh
git clone https://github.com/omyyys/diffuser.git
cd diffuser
make RAYLIB_DIR=./raylib/lib
./build/diffuser
```

On macOS, provide a compatible Raylib build in `raylib/lib`, or set `RAYLIB_DIR` to its library directory. The Makefile uses Clang-specific optimization flags, so a compatible toolchain is required.

## Usage

Select **Start** from the main menu. Hold the **left mouse button** and drag across the canvas to introduce fluid and movement.

| Control | Action |
| --- | --- |
| Left mouse drag | Add fluid and velocity |
| Left / Right arrow | Switch between density, velocity and vorticity views |
| `R` | Cycle through grid resolutions |
| `C` | Clear the simulation |
| `Esc` | Pause or resume |

The panel on the right lets you adjust **viscosity** and **diffusion**, choose a **colour**, and toggle **additive blending**. **Auto Colours** changes the colour over time, while **Auto Spawn** periodically introduces fluid at random positions; its interval can also be adjusted.

Close the application window to exit the simulation.

## Screenshots

![diffuser screenshot](https://github.com/user-attachments/assets/c8161239-3558-47af-bd39-bb7395d42b75)

![diffuser screenshot](https://github.com/user-attachments/assets/fc3e8a72-ea2b-4249-99c5-c72de69c1cd8)
