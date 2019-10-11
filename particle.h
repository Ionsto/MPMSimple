#pragma once
#include "glm/glm.hpp"
struct Particle{
    float Mass = 1;
    float E = 1e4;
    float v = 0;
    glm::vec2 Position;
    glm::vec2 Acceleration;
    glm::vec2 Velocity;
    glm::vec2 Momentum;
    glm::vec2 Force;
    glm::mat2x2 StrainRate;
    glm::mat2x2 Strain;
    glm::mat2x2 Stress;
};
