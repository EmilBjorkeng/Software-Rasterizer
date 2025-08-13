#ifndef __TIME_UTILS_H__
#define __TIME_UTILS_H__

#include <SDL3_ttf/SDL_ttf.h>
#include <cstdint>

extern float DeltaTime;
extern float fps;
extern float MSPerFrame;

extern uint64_t PerfFrequency;
extern uint64_t FPSCounter;

void deltaTimeGetPerformance();
void deltaTime();
void drawFPS(SDL_Renderer *Renderer, TTF_Font *Font);

#endif