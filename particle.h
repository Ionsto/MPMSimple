#pragma once
#include "glm/glm.hpp"
struct Particle{
    float Mass = 1;
    float Volume = 1;
    glm::vec2 Position = glm::vec2(0);
    glm::vec3 Colour = glm::vec3(0);
    glm::vec2 Velocity = glm::vec2(0);
    glm::mat2x2 VelocityField = glm::mat2(0);
    glm::mat2x2 DeformationGradient = glm::mat2(1);
};
