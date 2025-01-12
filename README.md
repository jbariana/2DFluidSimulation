# Fluid Simulation with SDL

This project is a 2D particle-based fluid simulation created using SDL (Simple DirectMedia Layer) in C. The simulation features particles that move, interact with each other through collision and pressure, and respond to user input such as picking up and releasing particles with the mouse.

## Features
- **Particle System**: Simulates particles with positions, velocities, and accelerations.
- **Gravity**: Particles experience gravitational force that affects their movement.
- **Collisions**: Particles bounce off the screen boundaries and interact with each other through collisions.
- **Pressure Simulation**: A simple pressure system allows particles to push away from each other when they are too close.
- **Mouse Interaction**: The user can pick up particles by clicking and dragging them, simulating fluid manipulation.

## Requirements
- **SDL2**: This project uses SDL2 for graphics rendering. You can install SDL2 following the instructions from the official [SDL website](https://www.libsdl.org/).
- **C Compiler**: A C compiler (like GCC) is needed to build the project.

## Installation

1. **Clone the repository:**

2. **Install SDL2 (if not already installed):**

3. **Compile the code:**
   
4. **Run the simulation:**

## Usage

- **Picking up particles**: Left-click and hold the mouse to pick up particles. The particles will be attracted towards the mouse position.
- **Releasing particles**: Release the left-click to drop the particles. They will continue to interact with other particles and the environment.

## Code Structure

- **main.c**: Contains the main simulation loop, initialization, and event handling.
- **Particle Structure**: Represents each particle in the simulation with properties such as position, velocity, density, and pressure.
- **updateParticles**: Updates the physics of particles, including gravity, velocity, position, and interactions.
- **pickUpParticles**: Handles the particle pickup mechanism when the mouse is pressed.
- **drawParticle**: Renders each particle to the screen as a circle.

## Key Constants

- **PARTICLE_AMT**: Number of particles in the simulation.
- **PARTICLE_SIZE**: Radius of each particle.
- **SCREEN_X, SCREEN_Y**: Screen width and height in pixels.
- **GRAVITY**: The gravitational force applied to the particles.
- **TIME_STEP**: Time step for each simulation update (approx. 60 FPS).
- **BOUNCE_DAMPING**: Energy loss factor when particles bounce off the screen edges.
- **MAX_VELOCITY**: Maximum velocity a particle can have.
- **PICKUP_RADIUS**: Radius of the mouse-based pickup circle.
