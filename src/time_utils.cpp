#include "time_utils.h"
#include <SDL3/SDL.h>
#include <iostream>

float DeltaTime;
float fps;
float MSPerFrame;

uint64_t PerfFrequency;
uint64_t FPSCounter;

void deltaTimeGetPerformance() {
    PerfFrequency = SDL_GetPerformanceFrequency();
    FPSCounter = SDL_GetPerformanceCounter();
}

void deltaTime() {
    uint64_t EndCounter = SDL_GetPerformanceCounter();
    uint64_t ElapsedTime = EndCounter - FPSCounter;
    DeltaTime = (float)ElapsedTime / (float)PerfFrequency;

    fps = (float)PerfFrequency / (float)ElapsedTime;
    MSPerFrame = 1000.0f * DeltaTime;
}

void drawFPS(SDL_Renderer *Renderer, TTF_Font *Font) {
    char text[16];
    sprintf(text, "FPS: %d", (int)fps);

    SDL_Color color = {255, 255, 255, 255};
    int length = (int)strlen(text);
    SDL_Surface *surfaceMessage = TTF_RenderText_Solid(Font, text, length, color);
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

    SDL_FRect messageRect = {720.0f - (float)surfaceMessage->w - 5.0f, 5.0f, (float)surfaceMessage->w, (float)surfaceMessage->h};
    SDL_RenderTexture(Renderer, message, NULL, &messageRect);
    SDL_DestroySurface(surfaceMessage);
    SDL_DestroyTexture(message);
}