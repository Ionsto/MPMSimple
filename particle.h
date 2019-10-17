#pragma once
#include "glm/glm.hpp"
struct Particle{
    int Type = 0;
    float Mass = 1;
    float YoungsModulus = -1e5;
    float PoissonsRatio = 0.5;
    float Viscosity = 0;
    glm::vec2 Position = glm::vec2(0);
    glm::vec3 Colour = glm::vec3(0);
    glm::vec2 Acceleration = glm::vec2(0);
    glm::vec2 Velocity = glm::vec2(0);
    glm::vec2 Force = glm::vec2(0);
    glm::mat2x2 StrainRate = glm::mat2x2(0);
    glm::mat2x2 Strain = glm::mat2x2(0);
    glm::mat2x2 Stress = glm::mat2x2(0);
};
