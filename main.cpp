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
#include <sstream>
#include <fstream>
#include <filesystem>
#include "TinyPngOut.hpp"
namespace fs = std::filesystem;
//std::default_random_engine generator;
//std::uniform_real_distribution<float> distribution(-0.5,0.5);
constexpr int RealSize = 40;
constexpr static float GridDim = 0.2;
constexpr static float inertial_scalar_inv = 1.0/(0.25 * GridDim * GridDim);
constexpr static int GridSize = static_cast<int>(static_cast<float>(RealSize)/GridDim);

constexpr static int MaxTime = 300;
constexpr static float Delay = 1.0/60.0;
constexpr static float DeltaTime = Phase::DeltaTime;
static constexpr int SubSteps = static_cast<int>(Delay/DeltaTime);
static constexpr int Resolution = 20;
static constexpr int RenderSize = RealSize * Resolution;

glm::mat2 ModelElastic(Particle & particle)
{
    auto F = particle.DeformationGradient;
    auto J = glm::determinant(F);
    auto F_T = glm::transpose(F);
    auto F_inv_T = glm::inverse(F_T);
    auto F_minus_F_inv_T = F - F_inv_T; 

    auto P_term_0 = particle.elastic_mu * (F_minus_F_inv_T);
    auto P_term_1 = particle.elastic_lambda * std::log(J) * F_inv_T;
    auto P = P_term_0 + P_term_1;
    auto stress = (1.0f / J) * (P * F_T);
    auto mean = 0.5*(stress[0][0]+stress[1][1]);
    auto radii = std::sqrt(std::pow(0.5*(stress[0][0]+stress[1][1]),2) + std::pow(stress[0][1],2));
    float p1 = mean + radii;
    float p0 = mean - radii;
    float max_principal = 1e4;
    float max_tresca = 5e5;
    float elastic_limit = 1e4;
    float tresca = radii * 2;
    particle.Colour.r = 255;
    particle.Colour.g = 0;
//    particle.Colour.b = 255 * std::clamp(float(p1 / max_principal),0.0f,1.0f);
    particle.Colour.b = 255 * std::clamp(float(tresca / max_tresca),0.0f,1.0f);
    //if(p1 > max_principal){
    //particle.elastic_lambda = 11.5e5 * std::clamp(2 - (float(tresca/elastic_limit)),0.8f,1.0f);
    if(radii * 2 > max_tresca){
        particle.Colour.r = 0;
        particle.Colour.g = 255;
        particle.Colour.b = 0;
        stress[0][0] = std::min(0.0f,stress[0][0]);
        stress[1][1] = std::min(0.0f,stress[1][1]);
//        stress[0][1] = 0;
//        stress[1][0] = 0;
//        particle.DeformationGradient = glm::mat2x2(1);
    }
    return stress;
}
glm::mat2 ModelTresca(Particle & particle)
{
    auto F = particle.DeformationGradient;
    auto J = glm::determinant(F);
    auto F_T = glm::transpose(F);
    auto F_inv_T = glm::inverse(F_T);
    auto F_minus_F_inv_T = F - F_inv_T; 

    float elastic_mu = 1e4;
    float elastic_lambda = 2e6;
    float theta = 30;
    float adheasion = 2;

    auto P_term_0 = elastic_mu * (F_minus_F_inv_T);
    auto P_term_1 = elastic_lambda * std::log(J) * F_inv_T;
    auto P = P_term_0 + P_term_1;
    auto stress = (1.0f / J) * (P * F_T);
    auto mean = 0.5*(stress[0][0]+stress[1][1]);
    auto radii = std::sqrt(std::pow(0.5*(stress[0][0]+stress[1][1]),2) + std::pow(stress[0][1],2));
    stress[0][0] = std::min(adheasion,stress[0][0]);
    stress[1][1] = std::min(adheasion,stress[1][1]);

    float p3 = mean - radii;
    float p0 = mean + radii;

    float max_tresca = 5e4;
    float tresca = radii * 2;

    particle.Colour.r = 255;
    particle.Colour.g = 100;
    particle.Colour.b = 0;

    if(tresca > max_tresca || std::abs(stress[0][1]) > max_shear){
        particle.Colour.r = 0;
        particle.Colour.g = 255;
        particle.Colour.b = 0;
        stress[0][0] = std::min(0.0f,stress[0][0]);
        stress[1][1] = std::min(0.0f,stress[1][1]);
        stress[0][1] = std::copysign(std::min(max_shear,std::abs(stress[0][1])),stress[0][1]);
        stress[1][0] = std::copysign(std::min(max_shear,std::abs(stress[1][0])),stress[0][1]);
//        particle.DeformationGradient = glm::mat2x2(1);
    }
    return stress;
}
glm::mat2 ModelMohrColoumb(Particle & particle)
{
    auto F = particle.DeformationGradient;
    auto J = glm::determinant(F);
    auto F_T = glm::transpose(F);
    auto F_inv_T = glm::inverse(F_T);
    auto F_minus_F_inv_T = F - F_inv_T; 

    float elastic_mu = 1e4;
    float elastic_lambda = 2e6;
    float theta = 30;
    float adheasion = 2;

    auto P_term_0 = elastic_mu * (F_minus_F_inv_T);
    auto P_term_1 = elastic_lambda * std::log(J) * F_inv_T;
    auto P = P_term_0 + P_term_1;
    auto stress = (1.0f / J) * (P * F_T);
    auto mean = 0.5*(stress[0][0]+stress[1][1]);
    auto radii = std::sqrt(std::pow(0.5*(stress[0][0]+stress[1][1]),2) + std::pow(stress[0][1],2));
    stress[0][0] = std::min(adheasion,stress[0][0]);
    stress[1][1] = std::min(adheasion,stress[1][1]);

    float p0 = mean - radii;
    float max_principal = 1e4;
    float max_tresca = 5e4;
    float max_shear = 5e2;
    float elastic_limit = 1e3;
    float tresca = radii * 2;
    particle.Colour.r = 255;
    particle.Colour.g = 100;
    particle.Colour.b = 0;
    
    //float sigma = mean - (radii * std::sin(theta));
    //float tau = (radii * std::cos(theta));
    float tresca_shear = (mean * std::sin(theta)) + (adheasion*std::cos(theta));
//    particle.Colour.b = 255 * std::clamp(float(p1 / max_principal),0.0f,1.0f);
//    particle.Colour.b = 255 * std::clamp(float(tresca / max_tresca),0.0f,1.0f);
    //if(p1 > max_principal){
    //particle.elastic_lambda = 11.5e5 * std::clamp(2 - (float(tresca/elastic_limit)),0.8f,1.0f);
    if(radii * 2 > max_tresca || std::abs(stress[0][1]) > max_shear){
        particle.Colour.r = 0;
        particle.Colour.g = 255;
        particle.Colour.b = 0;
        stress[0][0] = std::min(0.0f,stress[0][0]);
        stress[1][1] = std::min(0.0f,stress[1][1]);
        stress[0][1] = std::copysign(std::min(max_shear,std::abs(stress[0][1])),stress[0][1]);
        stress[1][0] = std::copysign(std::min(max_shear,std::abs(stress[1][0])),stress[0][1]);
//        particle.DeformationGradient = glm::mat2x2(1);
    }
    return stress;
}
glm::mat2 ModelWater(Particle & particle){
    float eos_stiffness = 500;
    float eos_power = 4;
    float rest_density = 1;
    float density = GridDim * GridDim * (particle.Mass / particle.Volume);
    float pressure = std::max(-0.05f, eos_stiffness * (std::pow(density / rest_density, eos_power) - 1));
    glm::mat2x2 stress = glm::mat2x2(
                -pressure, 0, 
                    0, -pressure
            );
    particle.Colour.r = 0;
    particle.Colour.g = 0;
    particle.Colour.b = 255;// * std::clamp(float(rest_density / density),0.0f,1.0f);

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
glm::mat2 ModelAir(Particle & particle){
    float eos_stiffness = 400;
    float eos_power = 4;
    float rest_density = 0.01;
    float density = GridDim * GridDim * (particle.Mass / particle.Volume);
    float pressure = std::max(-0.1f, eos_stiffness * (std::pow(density / rest_density, eos_power) - 1));
    glm::mat2x2 stress = glm::mat2x2(
                -pressure, 0, 
                    0, -pressure
            );
    particle.Colour.r = 0;
    particle.Colour.g = 255;
    particle.Colour.b = 0;

    float dynamic_viscosity = 0.0005;
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
Phase PhaseSoil;
Phase PhaseAir;

std::vector<uint8_t> frame(RenderSize * RenderSize * 4, 0);
std::vector<uint8_t> frame_png(RenderSize * RenderSize * 3, 0);

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
    size /= 2;
    for(int dx = -size;dx <= size;++dx){
        for(int dy = -size;dy <= size;++dy){
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
void SaveParticles(Phase & ph,GifWriter & g,int delay,int size = 3){
        for(int i = 0;i < ph.ParticleList.ParticleCount;++i){
            auto p = ph.ParticleList.Get(i);
            PaintXY(p.Position.x,RealSize - p.Position.y,p.Colour.r,p.Colour.g,p.Colour.b,size);
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
            //PaintVector(start,g.Velocity);
        }
    }
}
void PhaseCouplingNoSlip(Phase & PhA,Phase & PhB)
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
void PhaseCouplingFluid(Phase & PhA,Phase & PhB)
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
                auto drag = dv * 0.2f;
                GA.Velocity -= (drag * nB);
                GB.Velocity += (drag * nA);
            }
        }
    }
}
int main(int argc, char ** args)
{
	auto FramePath =fs::current_path() / "Frames/"; 
    fs::remove_all(FramePath);
    fs::create_directories(FramePath);
	std::cout<<FramePath<<std::endl;
	PhaseWater.model = ModelWater;
    PhaseElastic.model = ModelElastic;
    PhaseSoil.model = ModelMohrColoumb;
    PhaseAir.model = ModelAir;
	auto fileName = "out.gif";
	int delay = Delay;//round((DeltaTime * float(SubSteps))/1e-2);
    std::cout<<"Delay: "<<delay<<"\n";
    float density = 1.5;
    GifWriter g;
	GifBegin(&g, fileName, RenderSize, RenderSize, delay);
    using namespace std::chrono;
    auto start = high_resolution_clock::now();
    float WaterHeight = RealSize/5;
    //Setup
//void CreateRect(glm::vec2 pos,glm::vec2 size,float density = 40,float resolution = 0.1*1.5){
    //PhaseAir.CreateRect(glm::vec2(RealSize/2,RealSize/2),glm::vec2(RealSize/2,RealSize/2),20,0.3);
    //PhaseAir.CreateRectFixedMass(glm::vec2(RealSize/2,RealSize/2),glm::vec2(RealSize/2,RealSize/2),0.2,0.01);

    //PhaseWater.CreatePond(WaterHeight,0.18,1);
    //PhaseElastic.CreateBoat(glm::vec2(30,WaterHeight + 3));
    //PhaseElastic.CreateBoat(glm::vec2(15,WaterHeight + 3));
    float BeamLength = RealSize/3;
    //PhaseAir.CreateRectFixedMass(glm::vec2(RealSize/2,RealSize/2),glm::vec2(RealSize/2.0,RealSize/2.0),0.5,0.1);
    
    //PhaseElastic.CreateRect(glm::vec2(RealSize/2,8),glm::vec2(BeamLength + 2,0.5),40);
    //PhaseElastic.CreateRect(glm::vec2(RealSize/2 - BeamLength,4),glm::vec2(1,4));
    //PhaseElastic.CreateRect(glm::vec2(RealSize/2 + BeamLength,4),glm::vec2(1,4));

   // PhaseElastic.CreateRect(glm::vec2(RealSize/2,30),glm::vec2(0.5,10),90);
   // PhaseElastic.CreateRect(glm::vec2(RealSize/2,20),glm::vec2(5,5),200);
   //
    PhaseSoil.CreateRect(glm::vec2(RealSize/2,2),glm::vec2(RealSize,2),50);

    PhaseSoil.CreateRect(glm::vec2(RealSize/2,5),glm::vec2(4,1),50);
    PhaseSoil.CreateRect(glm::vec2(RealSize/2,7),glm::vec2(3,1),50);
    PhaseSoil.CreateRect(glm::vec2(RealSize/2,9),glm::vec2(1,1),50);

    //PhaseWater.CreateRect(glm::vec2((RealSize/4)-2,6),glm::vec2((RealSize/4)-2,2),20,0.1,glm::vec2(0,0));

    //PhaseSoil.CreateRect(glm::vec2(RealSize/2,4),glm::vec2(4,4),10);
    for(int t = 0;t < MaxTime;++t)
    {
        float dt = DeltaTime;
        std::cout<<"T:"<<t<<std::endl;
        if(t == 100){
            //PhaseElastic.CreateRect(glm::vec2(RealSize/2,30),glm::vec2(1,4),90);
            //CreateBlock();
            //CreateBoat();
//            PhaseElastic.CreateBoat(glm::vec2(RealSize/4,30));
        }
		if(t % 1 == 0 && t < 10){
//			PhaseAir.CreateRectFixedMass(glm::vec2(RealSize/2,36),glm::vec2(RealSize/2.0,1),0.2,0.05);
		}
		if(t % 1 == 0){
			//PhaseWater.CreateRect(glm::vec2(2,15),glm::vec2(1,1),2,0.5,glm::vec2(0,-2));
        }
        for(int i = 0; i < SubSteps;++i){
//            PhaseElastic.UpdateBegin();
            PhaseWater.UpdateBegin();
            PhaseSoil.UpdateBegin();
//            PhaseAir.UpdateBegin();
//            PhaseCouplingNoSlip(PhaseWater,PhaseElastic);
//            PhaseCouplingNoSlip(PhaseSoil,PhaseElastic);
            PhaseCouplingFluid(PhaseSoil,PhaseWater);
//            PhaseCoupling(PhaseWater,PhaseAir);
//            PhaseCoupling(PhaseElastic,PhaseAir);
//            PhaseElastic.UpdateEnd();
            PhaseWater.UpdateEnd();
            PhaseSoil.UpdateEnd();
//            PhaseAir.UpdateEnd();
        }
        //std::cout<< "Timings\n";
        //std::cout<< "Reset grid:" << Time_ResetGrid<<"\n";
        //std::cout<< "P2G:" << Time_P2G<<"\n";
        //std::cout<< "Update nodes:" << Time_UpdateNodes<<"\n";
        //std::cout<< "G2P:" << Time_G2P<<"\n";
        //std::cout<< "Update particle:" << Time_UpdateParticles<<"\n";
        std::fill(frame.begin(),frame.end(),0);
        //PaintGrid(PhaseAir);
        SaveParticles(PhaseAir,g,delay,2);
        SaveParticles(PhaseWater,g,delay);
        SaveParticles(PhaseElastic,g,delay);
        SaveParticles(PhaseSoil,g,delay);
        PaintNumbers(MaxTime,t);
        GifWriteFrame(&g, frame.data(), RenderSize,RenderSize, delay);
		std::ostringstream stringStream;
		stringStream << "frame " << t << " .png";
		auto filename = stringStream.str();
		auto path = FramePath/filename;

		for(int i = 0,j = 0;i < frame.size();){
			frame_png[j++] = frame[i++];
			frame_png[j++] = frame[i++];
			frame_png[j++] = frame[i++];
			++i;
		}
		std::ofstream out(path, std::ios::binary);
		TinyPngOut pngout(static_cast<uint32_t>(RenderSize), static_cast<uint32_t>(RenderSize), out);
		pngout.write(frame_png.data(), static_cast<size_t>(RenderSize * RenderSize));
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::fill(frame.begin(),frame.end(),255);
    GifWriteFrame(&g, frame.data(), RenderSize, RenderSize, delay); 
	GifEnd(&g);

    std::cout<<"Finished"<<std::endl;
    auto dur = duration_cast<duration<double>>(end-start);
    std::cout<<"Took "<< dur.count() <<"s"<<std::endl;
}
