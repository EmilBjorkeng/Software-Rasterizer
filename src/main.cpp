#define _USE_MATH_DEFINES
#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <iostream>
#include <string.h>
#include <cstdint>

#include "vector.h"
#include "matrix.h"
#include "matrix_vector_bridge.h"
#include "render.h"
#include "obj.h"
#include "camera.h"
#include "input.h"
#include "time_utils.h"
#include "shader.h"
#include "light.h"

#define ASSERT(_e, ...) if (!(_e)) { fprintf(stderr, __VA_ARGS__); exit(1); }

#define TITLE "Renderer"
#define WINDOW_WIDTH SCREEN_WIDTH
#define WINDOW_HEIGHT SCREEN_HEIGHT

int pitch;
uint32_t *drawBuffer;
float depthBuffer[WINDOW_WIDTH * WINDOW_HEIGHT];

SDL_Window *Window;
SDL_Renderer *Renderer;
SDL_Event WindowEvent;
TTF_Font *Font;
SDL_Texture *texture;

void render_text(char *text, int x, int y, SDL_Color color) {
    int i = 0;
    char *token = strtok(text, "\n");
    while (token != NULL)
    {
        int length = (int)strlen(token);
        SDL_Surface *surfaceMessage = TTF_RenderText_Solid(Font, token, length, color);
        if (!surfaceMessage) {
            SDL_Log("Failed to render text surface: %s", SDL_GetError());
            return;
        }
        SDL_Texture *message = SDL_CreateTextureFromSurface(Renderer, surfaceMessage);
        if (!message) {
            SDL_Log("Failed to create texture from surface: %s", SDL_GetError());
            SDL_DestroySurface(surfaceMessage);
            return;
        }
        SDL_FRect messageRect = { (float)x, (float)(y + i*40), (float)surfaceMessage->w, (float)surfaceMessage->h };
        SDL_RenderTexture(Renderer, message, NULL, &messageRect);
        SDL_DestroySurface(surfaceMessage);
        SDL_DestroyTexture(message);

        token = strtok(NULL, "\n");
        i++;
    }
    return;
}

int main(int argc, char *argv[]) {
    ASSERT(
        SDL_Init(SDL_INIT_EVENTS),
        "Initialization of SDL failed: %s\n",
        SDL_GetError());

    ASSERT(
        TTF_Init(),
        "Initialization of TTF failed: %s\n",
        SDL_GetError());

    Window = SDL_CreateWindow(
        TITLE,
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        SDL_WINDOW_HIGH_PIXEL_DENSITY);
    ASSERT(Window, "Window creation failed: %s\n", SDL_GetError());

    Renderer = SDL_CreateRenderer(Window, NULL);
    ASSERT(Renderer, "Creation of the renderer failed: %s\n", SDL_GetError());

    Font = TTF_OpenFont(
        "assets/Helvetica.ttf",
        24);
    ASSERT(Font, "Cound not load the font: %s\n", SDL_GetError());

    SDL_Texture* drawTexture = SDL_CreateTexture(Renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, WINDOW_WIDTH, WINDOW_HEIGHT);

    Camera<float> camera(float3(0, 2, 0), float3(0, 0, 0));

    std::vector<const Object*> sceneObjects;
    std::vector<const Light*> sceneLight;

    Object WorldAxis = LoadObjFile("assets/WorldAxis.obj");
    ASSERT(WorldAxis.faces.size() > 0, "Could not load assets/WorldAxis.obj\n");
    WorldAxis.scale = float3::Filled(0.2f);
    sceneObjects.push_back(&WorldAxis);

    Object CapriCube = LoadObjFile("assets/Cube.obj");
    CapriCube.position = {0, -1, 5};
    CapriCube.scale = float3::Filled(0.5f);
    ASSERT(CapriCube.faces.size() > 0, "Could not load assets/Cube.obj\n");
    CapriCube.shader = std::make_unique<TextureShader>();
    sceneObjects.push_back(&CapriCube);

    Object AlphaCube = LoadObjFile("assets/AlphaCube.obj");
    AlphaCube.position = {-2, 0, 4};
    AlphaCube.scale = float3::Filled(0.5f);
    ASSERT(AlphaCube.faces.size() > 0, "Could not load assets/AlphaCube.obj\n");
    AlphaCube.shader = std::make_unique<DiffuseAlphaShader>();
    sceneObjects.push_back(&AlphaCube);

    Object Monkey = LoadObjFile("assets/Monkey.obj");
    Monkey.position = {5, 0, 7};
    ASSERT(Monkey.faces.size() > 0, "Could not load assets/Monkey.obj\n");
    Monkey.shader = std::make_unique<LitDiffuseShader>();
    sceneObjects.push_back(&Monkey);

    Object LightObj = LoadObjFile("assets/SmallCube.obj");
    ASSERT(LightObj.faces.size() > 0, "Could not load assets/SmallCube.obj\n");
    LightObj.position = {5, -1, 6};
    LightObj.scale = float3::Filled(0.15f);
    LightObj.shader = std::make_unique<DiffuseShader>();
    sceneObjects.push_back(&LightObj);

    PointLight pointLight(float3(5, -1, 6), 1.0f, float3(0.0f, 1.0f, 0.0f));
    sceneLight.push_back(&pointLight);

    while (1) {
        deltaTimeGetPerformance();

        if (SDL_PollEvent(&WindowEvent)) {
            if (WindowEvent.type == SDL_EVENT_QUIT)
            { break; }
        }

        //
        // Input
        //
        const bool *keyboard_state_array = SDL_GetKeyboardState(NULL);

        // Camera Movement
        float moveSpeed = 2.0f;

        // Move forward/backward
        if (KEY_PRESSED(KEY_W)) {
            camera.position += camera.Forward() * moveSpeed * DeltaTime;
        }
        if (KEY_PRESSED(KEY_S)) {
            camera.position -= camera.Forward() * moveSpeed * DeltaTime;
        }
        // Move left/right
        if (KEY_PRESSED(KEY_A)) {
            camera.position -= camera.Right() * moveSpeed * DeltaTime;
        }
        if (KEY_PRESSED(KEY_D)) {
            camera.position += camera.Right() * moveSpeed * DeltaTime;
        }

        // Move up/down (world up)
        if (KEY_PRESSED(KEY_SPACE)) {
            camera.position.y() += moveSpeed * DeltaTime;
        }
        if (KEY_PRESSED(KEY_LCTRL) || KEY_PRESSED(KEY_LSHIFT)) {
            camera.position.y() -= moveSpeed * DeltaTime;
        }

        // Camera Rotation
        float turnSpeed = 2.0f;
        float angle = turnSpeed * DeltaTime;

        if (KEY_PRESSED(KEY_UP)) {
            float3 rightAxis = camera.Right();
            Quaternion<float> pitchRot = Quaternion<float>::FromAxisAngle(rightAxis, angle);
            camera.rotation = pitchRot * camera.rotation;
        }
        if (KEY_PRESSED(KEY_DOWN)) {
            float3 rightAxis = camera.Right();
            Quaternion<float> pitchRot = Quaternion<float>::FromAxisAngle(rightAxis, -angle);
            camera.rotation = pitchRot * camera.rotation;
        }
        if (KEY_PRESSED(KEY_LEFT)) {
            float3 worldUp(0, 1, 0);
            Quaternion<float> yawRot = Quaternion<float>::FromAxisAngle(worldUp, -angle);
            camera.rotation = yawRot * camera.rotation;
        }
        if (KEY_PRESSED(KEY_RIGHT)) {
            float3 worldUp(0, 1, 0);
            Quaternion<float> yawRot = Quaternion<float>::FromAxisAngle(worldUp, angle);
            camera.rotation = yawRot * camera.rotation;
        }
        if (KEY_PRESSED(KEY_Q)) {
            float3 forwardAxis = camera.Forward();
            Quaternion<float> rollRot = Quaternion<float>::FromAxisAngle(forwardAxis, -angle);
            camera.rotation = rollRot * camera.rotation;
        }
        if (KEY_PRESSED(KEY_E)) {
            float3 forwardAxis = camera.Forward();
            Quaternion<float> rollRot = Quaternion<float>::FromAxisAngle(forwardAxis, angle);
            camera.rotation = rollRot * camera.rotation;
        }

        // Normalize rotation after updates to prevent drift
        camera.rotation = camera.rotation.normalized();

        //Monkey.rotation.y() += DeltaTime;

        //
        // Logic
        //

        //
        // Draw
        //

        SDL_LockTexture(drawTexture, NULL, (void**)&drawBuffer, &pitch);

        // Clear screen
        for (int y = 0; y < WINDOW_HEIGHT; y++) {
            for (int x = 0; x < WINDOW_WIDTH; x++) {
                drawBuffer[y * (pitch / sizeof(uint32_t)) + x] = 0xFF;
                depthBuffer[(y * WINDOW_WIDTH) + x] = INFINITY;
            }
        }

        RendersceneObjects(drawBuffer, pitch, depthBuffer, camera, sceneObjects, sceneLight);

        /*renderObj(drawBuffer, depthBuffer, camera, WorldAxis);
        //renderObj(drawBuffer, depthBuffer, camera, Floor);
        renderObj(drawBuffer, depthBuffer, camera, Cube);
        renderObj(drawBuffer, depthBuffer, camera, Monkey);
        renderObj(drawBuffer, depthBuffer, camera, Dog);*/

        SDL_UnlockTexture(drawTexture);
        //SDL_UpdateTexture(drawTexture, NULL, drawBuffer, WINDOW_WIDTH * sizeof(Uint32));
        SDL_RenderTexture(Renderer, drawTexture, NULL, NULL);

        drawFPS(Renderer, Font);

        // Debug
        /*char text[128];
        float3 camRotation = camera.rotation.ToEulerAngles();
        sprintf(text, "--Camera--\nX: %f\nY: %f\nZ: %f\nPitch: %f\nYaw: %f\nRoll: %f", camera.position[0], camera.position[1], camera.position[2], camRotation[0], camRotation[1], camRotation[2]);
        render_text(text, 5, 5, (SDL_Color){255, 255, 255});*/

        SDL_RenderPresent(Renderer);
        deltaTime();
    }
    SDL_DestroyTexture(drawTexture);
    SDL_DestroyWindow(Window);
    SDL_DestroyRenderer(Renderer);
    return 0;
}