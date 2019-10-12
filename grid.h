#include "glm/glm.hpp"
struct Grid{
    float Mass;
    glm::vec2 Velocity = glm::vec2(0);
    glm::vec2 Force = glm::vec2(0);
    glm::vec2 Acceleration = glm::vec2(0);
};
