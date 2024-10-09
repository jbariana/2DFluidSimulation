#include <SDL.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


// Declaring constants
#define PARTICLE_AMT 1500 // Amount of particles
#define PARTICLE_SIZE 6   // Radius of particles
#define SCREEN_X 600       // Screen size x
#define SCREEN_Y 400       // Screen size y
#define GRAVITY 25.0f     // Gravitational acceleration (m/s²)
#define TIME_STEP 0.16f   // Time step for updates (approx. 60 FPS)
#define BOUNCE_DAMPING 0.5f // Energy loss factor on bounce
#define MAX_VELOCITY 200.0f // Maximum particle velocity
#define PRESSURE_COEFFICIENT 1.0f 
#define PICKUP_RADIUS 50.0f // Radius of the pickup circle


// Particle structure
typedef struct {
    float x, y; // Position
    float vx, vy; // Velocity
    float ax, ay; // Acceleration
    float density;
    float pressure;
    float radius; // Added radius
    int isClumped; // Indicates if the particle is part of a clump
} Particle;


// Declaring functions
Particle* populateParticles(int count);
Particle createParticle(float x, float y, float vx, float vy, float ax, float ay, float density, float pressure, float radius);
void drawParticle(SDL_Renderer* renderer, Particle particle);
void updateParticles(Particle* particles, int count);
void pickUpParticles(Particle* particles, int count, float mouseX, float mouseY);
void releaseParticles(Particle* particles, int count);
void applyMouseMovement(Particle* particles, int count, float mouseX, float mouseY);
void drawPickupCircle(SDL_Renderer* renderer, float mouseX, float mouseY);

// Main function
int main() {
    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return 1;
    }

    // Create an SDL window
    SDL_Window* window = SDL_CreateWindow("Fluid Simulator", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCREEN_X, SCREEN_Y, SDL_WINDOW_SHOWN);

    // Create a renderer
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    // Set the background color
    SDL_SetRenderDrawColor(renderer, 72, 0, 92, 255);

    // Clear the screen
    SDL_RenderClear(renderer);

    // Populate array of particles
    Particle* particles = populateParticles(PARTICLE_AMT);

    // Main loop
    int running = 1;
    int isPicking = 0; // Flag to check if particles are being picked up
    while (running) {
        // Clear the screen
        SDL_SetRenderDrawColor(renderer, 72, 0, 92, 255); // Background color
        SDL_RenderClear(renderer);

        // Get mouse state
        int mouseX, mouseY;
        Uint32 mouseState = SDL_GetMouseState(&mouseX, &mouseY);

        // Draw the pickup circle
        drawPickupCircle(renderer, mouseX, mouseY);

        //if mouse is clicked
        if (mouseState & SDL_BUTTON(SDL_BUTTON_LEFT)) {
            if (!isPicking) {
                // Start picking up particles
                pickUpParticles(particles, PARTICLE_AMT, mouseX, mouseY);
                isPicking = 1;
            }
            else {
                // While picking up, apply mouse movement
                applyMouseMovement(particles, PARTICLE_AMT, mouseX, mouseY);
            }
        }
        else {
            // Release particles when mouse is not pressed
            if (isPicking) {
                releaseParticles(particles, PARTICLE_AMT);
                isPicking = 0;
            }
        }

        // Update particles' physics
        updateParticles(particles, PARTICLE_AMT);

        // Draw particles
        for (int i = 0; i < PARTICLE_AMT; i++) {
            drawParticle(renderer, particles[i]);
        }

        // Draw the water surface

        // Present the renderer
        SDL_RenderPresent(renderer);

        // Closes if key is pressed
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            // Check for quit event (close button)
            if (event.type == SDL_QUIT) {
                running = 0;
            }
            // Check for key press event
            if (event.type == SDL_KEYDOWN) {
                running = 0;
            }
        }
    }

    // Frees allocated memory for particles
    free(particles);

    // Cleanup and quit SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}


// Draws particles
void drawParticle(SDL_Renderer* renderer, Particle particle) {
    // Set render color for the particle
    SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);

    int centerX = (int)particle.x;
    int centerY = (int)particle.y;
    int radius = (int)particle.radius;

    // Draw the circle with additional padding
    for (int y = -radius - 2; y <= radius; y++) { // Increased range for more padding
        for (int x = -radius - 2; x <= radius; x++) { // Increased range for more padding
            if ((x * x + y * y) <= (radius) * (radius)) { // Increased buffer
                SDL_RenderDrawPoint(renderer, centerX + x, centerY + y);
            }
        }
    }
}


// Update particle physics with gravity and bounce
void updateParticles(Particle* particles, int count) {
    for (int i = 0; i < count; i++) {
        Particle* p = &particles[i];

        // Apply gravity to vertical velocity
        p->vy += GRAVITY * TIME_STEP;

        // Update position based on velocity
        p->x += p->vx * TIME_STEP;
        p->y += p->vy * TIME_STEP;

        // Limit velocity to avoid instability
        float speed = sqrt(p->vx * p->vx + p->vy * p->vy);
        if (speed > MAX_VELOCITY) {
            p->vx *= MAX_VELOCITY / speed;
            p->vy *= MAX_VELOCITY / speed;
        }

        // Apply viscosity damping
        p->vx *= 0.99f; // Damping factor for viscosity
        p->vy *= 0.99f;

        // Boundary conditions: horizontal
        if (p->x >= SCREEN_X - p->radius) {
            p->x = SCREEN_X - p->radius;
            p->vx = -p->vx * BOUNCE_DAMPING;
        }
        else if (p->x < p->radius) {
            p->x = p->radius;
            p->vx = -p->vx * BOUNCE_DAMPING;
        }

        // Boundary conditions: vertical
        if (p->y >= SCREEN_Y - p->radius) {
            p->y = SCREEN_Y - p->radius;
            p->vy = -p->vy * BOUNCE_DAMPING;
        }
        else if (p->y < p->radius) {
            p->y = p->radius;
            p->vy = -p->vy * BOUNCE_DAMPING;
        }
    }

    // Collision detection and response
    for (int iter = 0; iter < 5; iter++) {
        for (int i = 0; i < count; i++) {
            for (int j = i + 1; j < count; j++) {
                Particle* p1 = &particles[i];
                Particle* p2 = &particles[j];

                // Calculate the distance between particles
                float dx = p1->x - p2->x;
                float dy = p1->y - p2->y;
                float distance = sqrt(dx * dx + dy * dy);
                float minDistance = p1->radius + p2->radius;

                // Check for collision
                if (distance < minDistance) {
                    // Collision response: move particles apart
                    float overlap = minDistance - distance;

                    // Normalized direction vector between particles
                    float nx = dx / distance;
                    float ny = dy / distance;

                    // Apply separation based on overlap
                    float separationFactor = overlap * 0.5f;

                    p1->x += nx * separationFactor;
                    p1->y += ny * separationFactor;
                    p2->x -= nx * separationFactor;
                    p2->y -= ny * separationFactor;

                    // Calculate the relative velocity
                    float relativeVx = p1->vx - p2->vx;
                    float relativeVy = p1->vy - p2->vy;

                    // Apply a basic impulse to simulate fluid pressure
                    float impulse = (relativeVx * nx + relativeVy * ny) * 0.5f; // Lower value for softer collisions

                    // Apply impulse to both particles
                    p1->vx -= impulse * nx * 0.5f; // Share impulse based on mass
                    p1->vy -= impulse * ny * 0.5f;
                    p2->vx += impulse * nx * 0.5f;
                    p2->vy += impulse * ny * 0.5f;
                }
            }
        }
    }

    for (int i = 0; i < count; i++) {
        Particle* p1 = &particles[i];
        p1->pressure = 0.0f; // Reset pressure

        for (int j = 0; j < count; j++) {
            if (i != j) {
                Particle* p2 = &particles[j];

                float dx = p2->x - p1->x;
                float dy = p2->y - p1->y;
                float distance = sqrt(dx * dx + dy * dy);
                float minDistance = p1->radius + p2->radius;

                if (distance < minDistance) {
                    float overlap = minDistance - distance;
                    p1->pressure += overlap; // Accumulate pressure
                }
            }
        }

        // Apply pressure as a repulsive force
        p1->vx += (p1->pressure / 100.0f) * (p1->vx < 0 ? 1 : -1); // Simple pressure response
        p1->vy += (p1->pressure / 100.0f) * (p1->vy < 0 ? 1 : -1);
    }

}


// Populates particles
Particle* populateParticles(int count) {
    Particle* particles = (Particle*)malloc(count * sizeof(Particle));
    for (int i = 0; i < count; i++) {
        float x = rand() % (SCREEN_X - 2 * PARTICLE_SIZE) + PARTICLE_SIZE;
        float y = rand() % (SCREEN_Y - 2 * PARTICLE_SIZE) + PARTICLE_SIZE;
        float vx = ((float)(rand() % 100) / 50.0f) - 1.0f; // Random velocity between -1 and 1
        float vy = ((float)(rand() % 100) / 50.0f) - 1.0f; // Random velocity between -1 and 1
        float ax = 0.0f; // Reset acceleration
        float ay = 0.0f; // Reset acceleration
        float density = 1.0f; // Default density
        float pressure = 0.0f; // Default pressure
        float radius = PARTICLE_SIZE; // Use constant radius
        particles[i] = createParticle(x, y, vx, vy, ax, ay, density, pressure, radius);
    }
    return particles;
}


// Creates a particle
Particle createParticle(float x, float y, float vx, float vy, float ax, float ay, float density, float pressure, float radius) {
    Particle p;
    p.x = x;
    p.y = y;
    p.vx = vx;
    p.vy = vy;
    p.ax = ax;
    p.ay = ay;
    p.density = density;
    p.pressure = pressure;
    p.radius = radius; // Store radius in the particle
    p.isClumped = 0; // Initially not clumped
    return p;
}


// Picks up particles when the mouse is pressed
void pickUpParticles(Particle* particles, int count, float mouseX, float mouseY) {
    const float PICKUP_RADIUS_SQUARED = PICKUP_RADIUS * PICKUP_RADIUS; // Pre-compute squared radius

    for (int i = 0; i < count; i++) {
        float dx = particles[i].x - mouseX;
        float dy = particles[i].y - mouseY;
        float distanceSquared = dx * dx + dy * dy; // Calculate squared distance

        if (distanceSquared < PICKUP_RADIUS_SQUARED) {
            particles[i].isClumped = 1; // Mark particle as clumped

            // Calculate the distance
            float distance = sqrt(distanceSquared); // Get the actual distance
            if (distance > 0) {
                float force = 0.8f; // Adjust this value for sensitivity
                float scale = force / distance; // Normalize the force by distance

                // Add velocity towards the mouse position
                particles[i].vx += dx * scale; // Smoothly move towards mouse
                particles[i].vy += dy * scale; // Smoothly move towards mouse
            }
        }
    }
}


// Releases particles when the mouse is not pressed
void releaseParticles(Particle* particles, int count) {
    for (int i = 0; i < count; i++) {
        particles[i].isClumped = 0; // Mark particle as not clumped
    }
}


// Applies mouse movement to clumped particles
void applyMouseMovement(Particle* particles, int count, float mouseX, float mouseY) {
    const float separationDistance = 2.0f * PARTICLE_SIZE; // Minimum allowed distance between clumped particles
    const float separationForceStrength = 30.0f; // Strength of the separation force
    const float pickupRadius = 50.0f; // Radius within which particles will not move towards the mouse

    for (int i = 0; i < count; i++) {
        if (particles[i].isClumped) {
            // Calculate distance to mouse
            float dx = mouseX - particles[i].x;
            float dy = mouseY - particles[i].y;
            float distance = sqrtf(dx * dx + dy * dy);

            // Only apply movement if the particle is outside the pickup radius
            if (distance > pickupRadius) {
                float force = 30.0f; // Adjust this value for sensitivity
                float scale = force / distance; // Normalize the force by distance

                // Smoothly move towards mouse
                particles[i].vx += (dx * scale);
                particles[i].vy += (dy * scale);

                // Decay to slow down clumped particles
                particles[i].vx *= 0.96f;
                particles[i].vy *= 0.96f;
            }

            // Apply separation force to keep clumped particles from overlapping
            for (int j = 0; j < count; j++) {
                if (i != j && particles[j].isClumped) {
                    float dx = particles[i].x - particles[j].x;
                    float dy = particles[i].y - particles[j].y;
                    float distance = sqrtf(dx * dx + dy * dy);

                    // If the particles are too close, apply a separation force
                    if (distance < separationDistance && distance > 0) {
                        float separationForce = separationForceStrength / distance; // Stronger when closer
                        float nx = dx / distance; // Normalize the vector
                        float ny = dy / distance;

                        // Push the particles away from each other
                        particles[i].vx += nx * separationForce;
                        particles[i].vy += ny * separationForce;
                        particles[j].vx -= nx * separationForce;
                        particles[j].vy -= ny * separationForce;
                    }
                }
            }
        }
    }
}


// Draws the pickup circle at the mouse position
void drawPickupCircle(SDL_Renderer* renderer, float mouseX, float mouseY) {
    // Set the color for the pickup circle
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

    // Draw a circle by using lines
    for (float angle = 0; angle < 360; angle += 0.1f) {
        float radian = angle * (M_PI / 180); // Convert degrees to radians
        int x = (int)(mouseX + PICKUP_RADIUS * cos(radian)); // Calculate x
        int y = (int)(mouseY + PICKUP_RADIUS * sin(radian)); // Calculate y
        SDL_RenderDrawPoint(renderer, x, y); // Draw the point on the circle
    }
}

__declspec(dllexport) void startSimulation() {
    main();
}