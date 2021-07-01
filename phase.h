#pragma once
#include <iostream>
#include "particle.h"
#include "grid.h"
#include "gif.h"
#include "SwapList.h"
#include <chrono>
#include <array>
#include <bitset>
#include <vector>
#include <random>
#include <math.h>
#include <algorithm>
struct Phase{
    std::function<glm::mat2(Particle &)> model;
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(-0.5,0.5);
    constexpr static float DeltaTime = 5e-4;
    constexpr static int MaxParticles = 100000;
    constexpr static int RealSize = 40;
    constexpr static float GridDim = 0.2;
    constexpr static float inertial_scalar_inv = 1.0/(0.25 * GridDim * GridDim);
    constexpr static int GridSize = static_cast<int>(static_cast<float>(RealSize)/GridDim);
    SwapList<Particle,MaxParticles> ParticleList;
    std::array<Grid,GridSize * GridSize> SimGrid;
    Grid & GetGrid(int x, int y){
        return SimGrid[(x * GridSize) + y];
    }
    float WeightAxis(float x)
    {
        float eta = x/GridDim;
        if((eta < -0.5) && (eta > -1.5))
        {
            return 0.125 * std::pow(3 + (2 * eta),2);
        }
        if((eta > 0.5) && (eta < 1.5))
        {
            return 0.125 * std::pow(3 - (2 * eta),2);
        }
        if((eta >= -0.5) && (eta <= 0.5))
        {
            return 0.75 - std::pow(eta,2);
        }
        return 0;
    }
    float Weight(glm::vec2 distance)
    {
        return WeightAxis(distance.x) * WeightAxis(distance.y);
    }
    template <typename F>
    void IterateOverNeighbours(Particle & p,F f)
    {
        int GridX = std::round(p.Position.x / GridDim) - 1;
        int GridY = std::round(p.Position.y / GridDim) - 1;
        for(int dx = 0;dx < 3;++dx)
        {
            for(int dy = 0;dy < 3;++dy)
            {
                int x = GridX + dx;
                int y = GridY + dy;
                if(x >= 0 && x < GridSize){
                    if(y >= 0 && y < GridSize){
                        auto d = p.Position - (glm::vec2(static_cast<float>(x),static_cast<float>(y)) * GridDim);
                        float weight = Weight(d);
                        f(p,x,y); 
                    }
                }
            }
        }
    }
    glm::mat2 ModelElastic(Particle & particle)
    {
        auto F = particle.DeformationGradient;
        auto J = glm::determinant(F);
        auto F_T = glm::transpose(F);
        auto F_inv_T = glm::inverse(F_T);
        auto F_minus_F_inv_T = F - F_inv_T; 

        //Lame properties
        float elastic_mu = 8.2e5; 
        float elastic_lambda = 11.5e5;
        auto P_term_0 = elastic_mu * (F_minus_F_inv_T);
        auto P_term_1 = elastic_lambda * std::log(J) * F_inv_T;
        auto P = P_term_0 + P_term_1;
        auto stress = (1.0f / J) * (P * F_T);
        return stress;
    }
    glm::mat2 ModelWater(Particle & particle){
        float eos_stiffness = 10;
        float eos_power = 4;
        float rest_density = 0.3 * 2;
        float density = GridDim * GridDim * (particle.Mass / particle.Volume);
        float pressure = std::max(-0.1f, eos_stiffness * (std::pow(density / rest_density, eos_power) - 1));
        glm::mat2x2 stress = glm::mat2x2(
                    -pressure, 0, 
                        0, -pressure
                );
        particle.Colour.r = 0;
        particle.Colour.g = 255;
        particle.Colour.b = 255 * std::clamp(float(rest_density / density),0.0f,1.0f);

        float dynamic_viscosity = 0.01;
        // velocity gradient - MLS-MPM eq. 17, where derivative of quadratic polynomial is linear
        glm::mat2x2 dudv = particle.VelocityField;
         // build strain from the velocity gradient
        float trace = (dudv[0][0] + dudv[1][1]);
        glm::mat2x2 strain = dudv;
        strain[1][0] = trace;
        strain[0][1] = trace;
        glm::mat2x2 viscosity_term = dynamic_viscosity * strain;
        stress += viscosity_term;
        return stress;
    }
    void P2G()
    {
#pragma omp parallel for
        for(int i = 0; i < ParticleList.ParticleCount;++i)
        {
            auto & particle = ParticleList.Get(i);
            IterateOverNeighbours(particle,
                [] (auto & p, int x ,int y){
                //APIC
                Grid & g = GetGrid(x,y);
                auto d = (glm::vec2(static_cast<float>(x),static_cast<float>(y)) * GridDim) - p.Position;
                float weight = Weight(d);
                auto WeightedMass = p.Mass * weight;
                auto velocity_updated = WeightedMass * (p.Velocity + (p.VelocityField*d));
#pragma omp atomic update
                    g.Mass += WeightedMass;
#pragma omp atomic update
                    g.Velocity[0] += velocity_updated[0];
#pragma omp atomic update
                    g.Velocity[1] += velocity_updated[1];
                });
        }
        //Do constitutive model update
#pragma omp parallel for
        for(int i = 0; i < ParticleList.ParticleCount;++i)
        {
            auto & particle = ParticleList.Get(i);
            float density = 0;
            IterateOverNeighbours(particle,
                    [&density] (auto & p, int x ,int y){
                        Grid g = GetGrid(x,y);
                        auto d = (glm::vec2(static_cast<float>(x),static_cast<float>(y)) * GridDim) - p.Position;
                        float w = Weight(d);
                        density += g.Mass * w;
                    });
            particle.Volume = GridDim * GridDim * (particle.Mass / density);
            glm::mat2 stress = model(particle);

            auto eq_16_term_0 = -particle.Volume * inertial_scalar_inv * stress * DeltaTime;
            IterateOverNeighbours(particle,
                [eq_16_term_0] (auto & p, int x ,int y){
                    Grid & g = GetGrid(x,y);
                    auto d = (glm::vec2(static_cast<float>(x),static_cast<float>(y)) * GridDim) - p.Position;
                    float w = Weight(d);
                    glm::vec2 dv = w * eq_16_term_0 * d; 
                    for(int d = 0;d < 2;++d){
#pragma omp atomic update
                        g.Velocity[d] += dv[d];
                    }
                });
            }
    }
    void UpdateNodes()
    {
#pragma omp parallel for// schedule (dynamic)
        for(auto & g : SimGrid){
            if(g.Mass != 0)
            {
                g.Velocity /= g.Mass;
                auto Acceleration = glm::vec2(0,-9.8);
                g.Velocity += Acceleration * DeltaTime;
            }
        }
        float Friction = 0.01;
#pragma omp parallel for
        for(int x = 0; x < GridSize;++x){
            for(int d = 0;d < 3;++d){
                //Floor 
                GetGrid(x,d).Velocity.y = std::max(GetGrid(x,d).Velocity.y,0.0f);
                GetGrid(x,d).Velocity.x *= Friction;
                //Ceiling
                //GetGrid(x,GridSize - (1 + d)).Velocity.y = 0;
                GetGrid(x,GridSize - (1 + d)).Velocity.y = std::min(GetGrid(x,GridSize - (1 + d)).Velocity.y,0.0f);
                GetGrid(x,GridSize - (1 + d)).Velocity.x *= Friction;
                //Left wall
                //GetGrid(d,x).Velocity.x = 0; 
                GetGrid(d,x).Velocity.x = std::max(GetGrid(d,x).Velocity.x,0.0f); 
                GetGrid(d,x).Velocity.y *= Friction;
                //Right wall
                //GetGrid(GridSize - (1 + d),x).Velocity.x = 0;
                GetGrid(GridSize - (1 + d),x).Velocity.x = std::min(GetGrid(GridSize - (1 + d),x).Velocity.x,0.0f);
                GetGrid(GridSize - (1 + d),x).Velocity.y *= Friction;
            }
        }
    }
    void G2PNode(int x,int y,Particle & p){
        Grid g = GetGrid(x,y);
        auto d = (glm::vec2(static_cast<float>(x),static_cast<float>(y)) * GridDim) - p.Position;
        float weight = Weight(d);
        p.Velocity += g.Velocity * weight;
        p.VelocityField += inertial_scalar_inv * glm::outerProduct(g.Velocity,d) * weight;
    }

    void G2P()
    {
#pragma omp parallel for
        for(int i = 0; i < ParticleList.ParticleCount;++i)
        {
            auto & particle = ParticleList.Get(i);
            particle.Velocity = glm::vec2(0);
            particle.VelocityField = glm::mat2x2(0);
            IterateOverNeighbours(particle,
                    [] (auto & p, int x ,int y){
                        G2PNode(x,y,p);
                    });
            auto FpNew = glm::mat2(1);
            FpNew += DeltaTime * particle.VelocityField;
            particle.DeformationGradient = FpNew * particle.DeformationGradient;
           
        }
    }
    void UpdateParticles()
    {
#pragma omp parallel for
        for(int i = 0; i < ParticleList.ParticleCount;++i)
        {
            auto & particle = ParticleList.Get(i);
            particle.Position += particle.Velocity * DeltaTime;

            //const float K_water = 50;
            //const float Gamma_water = 3;
            //float dJp = -K_water * (1.0 / pow(particle.Jp, Gamma_water) - 1.0);
            //particle.Ap = dJp * particle.Volume * particle.Jp;
        }

        for(int i = 0; i < ParticleList.ParticleCount;++i)
        {
            auto & particle = ParticleList.Get(i);
            const float boarder = (GridDim * 3);
            particle.Position.x = std::clamp(particle.Position.x,boarder, RealSize - boarder);
            particle.Position.y = std::clamp(particle.Position.y,boarder, RealSize - boarder);
    //Check in bounds
            if(particle.Position.x < 0)
            {
                particle.Position.x = 0;
            }
            if(particle.Position.x >= RealSize)
            {
                particle.Position.x = RealSize - 0.1;
            }
            if(particle.Position.y < 0)
            {
                particle.Position.y = 0;
            }
            if(particle.Position.y >= RealSize)
            {
                particle.Position.y = RealSize - 0.1;
            }

            if(isnan(particle.Position.x) || isnan(particle.Position.y)){
                ParticleList.Remove(i);
                i -= 1;
                continue;
            }
        }
    }
    template<typename F>
    double Timeit(F f, std::string Name){
        using namespace std::chrono;
        auto start = high_resolution_clock::now();
        f();
        auto end = high_resolution_clock::now();
        auto dur = duration_cast<duration<double>>(end-start);
        return dur.count();
    }
    double Time_ResetGrid = 0;
    double Time_P2G = 0;
    double Time_UpdateNodes = 0;
    double Time_G2P = 0;
    double Time_UpdateParticles = 0;
    void Update()
    {
        Time_ResetGrid += Timeit([=](){
        std::fill(SimGrid.begin(),SimGrid.end(),Grid());
        },"Reset grid");
        Time_P2G += Timeit(P2G,"P2G");
        Time_UpdateNodes += Timeit(UpdateNodes,"Update nodes");
        Time_G2P += Timeit(G2P,"G2P");
        Time_UpdateParticles += Timeit(UpdateParticles,"Update particles");
    }
    void AddParticle(glm::vec2 pos){
        auto p = Particle();
        p.Position = pos;
        ParticleList.Add(p);
    }
void CreateRect(glm::vec2 pos,glm::vec2 size,float density = 40,float resolution = 0.1*1.5){
    float noise = 0.01f;
    glm::ivec2 count = 2.0f*size / resolution;
    float mass_per_unit = (density*(4 * size.x*size.y)) / (count.x * count.y);
    for(int x = 0;x <= count.x;++x)
    {
        for(int y = 0;y <= count.y;++y)
        {
            auto pa = Particle();
            float dx = -size.x + (x*2*size.x/count.x) + (noise * distribution(generator));
            float dy = -size.y + (y*2*size.y/count.y) + (noise * distribution(generator));
            pa.Position = pos + glm::vec2(dx,dy); 
            pa.Colour.r = 255;
            pa.Mass = mass_per_unit;
            pa.Type = 0;
            ParticleList.Add(pa);
        }
    }
}
void CreateBlock(){
    int count = 20;
    float size = 2;
    float noise = 0.1;
    for(int x = -count;x <= count;++x)
    {
    for(int y = -count;y <= count;++y)
    {
        auto pa = Particle();
        float dx = (x*size/count) + (noise * distribution(generator));
        float dy = (y*size/count) + (noise * distribution(generator));
//                pa.Position = (glm::vec2(3,RealSize / 2.0f)) + glm::vec2(dx,dy); 
        pa.Position = (glm::vec2(RealSize,RealSize) / 2.0f) + glm::vec2(dx,dy); 
        pa.Colour.r = rand()%255;
        pa.Colour.g = rand()%255;
        pa.Colour.b = rand()%255;
        pa.Mass = 0.2;
        pa.Type = 0;
        ParticleList.Add(pa);
    }
    }
}
void CreateBoat(glm::vec2 pos = (glm::vec2(RealSize,RealSize) / 2.0f)){
    CreateRect(pos                 ,glm::vec2(3,0.5));
    CreateRect(pos + glm::vec2(3.5,1.5),glm::vec2(0.5,2));
    CreateRect(pos + glm::vec2(-3.5,1.5),glm::vec2(0.5,2));
}
void CreatePond(float height,float resolution = 0.22)
{
    float rest_density = 0.3 * 2;
    float mass = 2;
    float noise = 0.01;
    const float border = (GridDim * 3);
    glm::vec2 count = glm::vec2(RealSize-(border*2),height-(border)) / resolution;
    for(int x = 0;x < count.x;++x)
    {
        for(int y = 0;y < count.y;++y)
        {
                auto pa = Particle();
                pa.Mass = 2;
                float dx = noise * distribution(generator);
                float dy = noise * distribution(generator);
                pa.Position = (glm::vec2(x,y)*resolution)  + glm::vec2(dx,dy) + glm::vec2(border,border); 
                pa.Velocity = 1.0f * glm::vec2(distribution(generator),distribution(generator));
                pa.Velocity.x += 5;
                pa.Colour.r = rand()%255;
                pa.Colour.g = rand()%255;
                pa.Colour.b = rand()%255;
                ParticleList.Add(pa);
        }
    }
}
float WaterFlowCount = 0;
int ParticleCount = 0;
void WaterFlow(glm::vec2 pos,glm::vec2 size,float flow){
        WaterFlowCount += flow * DeltaTime;
        for(int p = 0;p < WaterFlowCount;++p)
        {
            WaterFlowCount--;
            auto pa = Particle();
            pa.Mass = 2;
            float dx = size.x * distribution(generator);
            float dy = size.y * distribution(generator);
//                pa.Position = (glm::vec2(3,RealSize / 2.0f)) + glm::vec2(dx,dy); 
            pa.Position = pos + glm::vec2(dx,dy); 
            //pa.Velocity.x = (8*sinf(3.14 * (float(t) / 20.0)));
            pa.Velocity.y = -5;
            pa.Colour.r = rand()%255;
            pa.Colour.g = rand()%255;
            pa.Colour.b = rand()%255;
            ParticleList.Add(pa);
        }
}
}
