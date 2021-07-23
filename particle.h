#pragma once
#include "glm/glm.hpp"
struct Particle{
    int Type = 0;
    float Mass = 1;
    float Volume = 5;
    float elastic_mu = 8.2e5; 
    float elastic_lambda = 10.0e5;
    glm::vec2 Position = glm::vec2(0);
    glm::vec3 Colour = glm::vec3(0);
    glm::vec2 Velocity = glm::vec2(0);
    glm::mat2x2 VelocityField = glm::mat2(0);
    glm::mat2x2 DeformationGradient = glm::mat2(1);
};
