#include <iostream>
#include "particle.h"
#include "phase.h"
#include "grid.h"
#include "gif.h"
#include "SwapList.h"
#include <chrono>
#include <array>
#include <bitset>
#include <vector>
#include <random>
#include <algorithm>
#include <math.h>
//std::default_random_engine generator;
//std::uniform_real_distribution<float> distribution(-0.5,0.5);
constexpr int RealSize = 40;
constexpr static float GridDim = 0.2;
constexpr static float inertial_scalar_inv = 1.0/(0.25 * GridDim * GridDim);
constexpr static int GridSize = static_cast<int>(static_cast<float>(RealSize)/GridDim);

constexpr static int MaxTime = 200;
constexpr static int Delay = 2;
constexpr static float DeltaTime = 5e-4;
static constexpr int SubSteps = (Delay * 1e-2)/DeltaTime;
static constexpr int Resolution = 20;
static constexpr int RenderSize = RealSize * Resolution;

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

Phase PhaseWater;
Phase PhaseElastic;

std::vector<uint8_t> frame(RenderSize * RenderSize * 4, 0);

void PaintPixel(float xi,float yi,int r,int g,int b){
    int x = static_cast<int>((xi));
    int y = static_cast<int>((yi));
    int v = 4 * (x + (RenderSize * y));
    if(x >= 0 && x < RenderSize){
        if(y >= 0 && y < RenderSize){
            frame[v] = r;
            frame[v+1] = g;
            frame[v+2] = b;
            frame[v+3] = 255;
        }
    }
}
void PaintXY(float xi,float yi,int r,int g,int b,int size=Resolution){
    for(int dx = 0;dx < size;++dx){
        for(int dy = 0;dy < size;++dy){
            int x = static_cast<int>((xi) * Resolution) + dx;
            int y = static_cast<int>((yi) * Resolution) + dy;
            int v = 4 * (x + (RenderSize * y));
            if(x >= 0 && x < RenderSize){
                if(y >= 0 && y < RenderSize){
                    frame[v] = r;
                    frame[v+1] = g;
                    frame[v+2] = b;
                    frame[v+3] = 255;
                }
            }
        }
    }
}
void SaveParticles(Phase & ph,GifWriter & g,int delay){
        for(int i = 0;i < ph.ParticleList.ParticleCount;++i){
            auto p = ph.ParticleList.Get(i);
            PaintXY(p.Position.x,RealSize - p.Position.y,p.Colour.r,p.Colour.g,p.Colour.b,5);
        }
}
void PaintNumbers(int MaxTime,int t)
{
        float FillPC = (float)t / (float)MaxTime;
        static constexpr int Bits = 5;
        int fillint = static_cast<int>(FillPC * (float)(1<<(Bits-1)));
        std::bitset<Bits> pattern = std::bitset<Bits>(fillint);
        for(int b = 0;b < Bits;++b){
            if(pattern[b])
            {
                PaintXY(b,0,0,255,0);
            }
        } PaintXY(Bits-1,0,255,255,0);
}
void PaintVector(glm::vec2 start,glm::vec2 delta){
    start *= Resolution;
    delta *= Resolution;
    PaintPixel(start.x,start.y,0,0,255);
    if(delta.x != 0 || delta.y != 0){
        float steplength = 0.1;
        auto dvec = glm::normalize(delta) * steplength;
        int step = glm::length(delta) / steplength;
        for(int i = 1;i < step;++i)
        {
            auto pos = start + (float(i) * dvec);
            int x = pos.x;
            int y = RenderSize - pos.y;
            PaintPixel(x,y,0,0,255);
        }
    }
}
void PaintGrid(Phase & ph)
{
    for(int x = 0;x < GridSize;++x){
        for(int y = 0;y < RenderSize;++y){
            PaintPixel(x * GridDim * Resolution,y,255,0,0);
            PaintPixel(y,x * GridDim * Resolution,255,0,0);
        }
    }
    for(int x = 0;x < GridSize;++x){
        for(int y = 0;y < GridSize;++y){
            Grid & g = ph.GetGrid(x,y);
            auto start = (glm::vec2(static_cast<float>(x),static_cast<float>(y)) * GridDim);
            PaintVector(start,g.Velocity);
        }
    }
}
void PhaseCoupling(Phase & PhA,Phase & PhB)
{
    for(int x = 0; x < PhA.GridSize;++x)
    {
        for(int y = 0; y < PhA.GridSize;++y)
        {
            auto& GA = PhA.GetGrid(x,y);
            auto& GB = PhB.GetGrid(x,y);
            auto dv = GA.Velocity - GB.Velocity;
            float MassTotal = GA.Mass + GB.Mass;
            if(MassTotal != 0)
            {
                float nA = GA.Mass / MassTotal;
                float nB = GB.Mass / MassTotal;
                auto drag = dv * 1.0f;
                GA.Velocity -= (drag * nB);
                GB.Velocity += (drag * nA);
            }
        }
    }
}
int main(int argc, char ** args)
{

    PhaseWater.model = ModelWater;
    PhaseElastic.model = ModelElastic;
	auto fileName = "out.gif";
	int delay = round((DeltaTime * float(SubSteps))/1e-2);
    std::cout<<"Delay: "<<delay<<"\n";
    float density = 1.5;
    GifWriter g;
	GifBegin(&g, fileName, RenderSize, RenderSize, delay);
    using namespace std::chrono;
    auto start = high_resolution_clock::now();
    float WaterHeight = RealSize/5;
    PhaseWater.CreatePond(WaterHeight);
    PhaseElastic.CreateBoat(glm::vec2(10,WaterHeight + 3));
    PhaseElastic.CreateBoat(glm::vec2(30,WaterHeight + 3));
    for(int t = 0;t < MaxTime;++t)
    {
        float dt = DeltaTime;
        std::cout<<"T:"<<t<<std::endl;
        if(t == 120){
            //CreateBlock();
            //CreateBoat();
        }
        for(int i = 0; i < SubSteps;++i){

            if(t < 100){
//                WaterFlow(glm::vec2(6,3),glm::vec2(5,5),2000 * t / 20);
            }
            PhaseWater.UpdateBegin();
            PhaseElastic.UpdateBegin();
            PhaseCoupling(PhaseWater,PhaseElastic);
            PhaseWater.UpdateEnd();
            PhaseElastic.UpdateEnd();
        }
        //std::cout<< "Timings\n";
        //std::cout<< "Reset grid:" << Time_ResetGrid<<"\n";
        //std::cout<< "P2G:" << Time_P2G<<"\n";
        //std::cout<< "Update nodes:" << Time_UpdateNodes<<"\n";
        //std::cout<< "G2P:" << Time_G2P<<"\n";
        //std::cout<< "Update particle:" << Time_UpdateParticles<<"\n";
        std::fill(frame.begin(),frame.end(),0);
        //PaintGrid();
        SaveParticles(PhaseWater,g,delay);
        SaveParticles(PhaseElastic,g,delay);
        PaintNumbers(MaxTime,t);
        GifWriteFrame(&g, frame.data(), RenderSize,RenderSize, delay);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::fill(frame.begin(),frame.end(),255);
    GifWriteFrame(&g, frame.data(), RenderSize, RenderSize, delay); 
	GifEnd(&g);

    std::cout<<"Finished"<<std::endl;
    auto dur = duration_cast<duration<double>>(end-start);
    std::cout<<"Took "<< dur.count() <<"s"<<std::endl;
}
