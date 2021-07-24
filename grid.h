#pragma once
#include "glm/glm.hpp"
struct Grid{
    float Mass = 0;
    float Volume = 0;
    glm::vec2 Velocity = glm::vec2(0,0);
    glm::vec2 Force = glm::vec2(0,0);
};
