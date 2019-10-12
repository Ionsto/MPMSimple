#include <iostream>
#include "particle.h"
#include "grid.h"
#include "gif.h"
#include <array>
#include <bitset>
#include <vector>
constexpr static int MaxParticles = 1000;
constexpr static int GridSize = 100;
constexpr static float DeltaTime = 1e-3;
static constexpr int SubSteps = 100;
static constexpr int Resolution = 10;
static constexpr int RenderSize = GridSize * 10;
int ParticleCount = 0;
std::array<Particle,MaxParticles> ParticleList;
std::array<Grid,GridSize * GridSize> SimGrid;
std::vector<uint8_t> frame(RenderSize * RenderSize * 4, 0);
Grid & GetGrid(int x, int y){
    return SimGrid[(x * GridSize) + y];
}
float WeightAxis(float x)
{
    return 1.0 - std::abs(x);
}
float WeightGradAxis(float x)
{
    return std::copysign(-1,x);
}
float Weight(glm::vec2 distance)
{
    return WeightAxis(distance.x) * WeightAxis(distance.y);
}
glm::vec2 WeightGrad(glm::vec2 distance)
{
    return glm::vec2(WeightGradAxis(distance.x) * WeightAxis(distance.y),
    WeightGradAxis(distance.y) * WeightAxis(distance.x));
}
void APIC(int x,int y,Particle p){
    Grid & g = GetGrid(x,y);
    auto d = p.Position - glm::vec2(static_cast<float>(x),static_cast<float>(y));
    float weight = Weight(d);
    auto weightgrad = WeightGrad(d);
    g.Mass += p.Mass * weight;
    g.Velocity += p.Velocity * p.Mass * weight;
    g.Force += (p.Stress * weightgrad);
    g.Force += p.Force * weight;
}
void P2G()
{
    for(int i = 0; i < ParticleCount;++i)
    {
        auto particle = ParticleList[i];
        int GridX = std::floor(particle.Position.x);
        int GridY = std::floor(particle.Position.y);
        for(int x = 0;x < 2;++x)
        {
            for(int y = 0;y < 2;++y)
            {
                APIC(GridX+x,GridY+y,particle);
            }
        }
    }
    for(auto & g : SimGrid){
        if(g.Mass != 0)
        {
            g.Velocity /= g.Mass;
            g.Acceleration = g.Force / g.Mass;
        }
    }
    for(int x = 0; x < GridSize;++x){
        GetGrid(x,1).Velocity.y = 0;
        GetGrid(x,0).Velocity.y = 0;
    }
}
void G2PNode(int x,int y,Particle & p){
    Grid g = GetGrid(x,y);
    auto d = p.Position - glm::vec2(static_cast<float>(x),static_cast<float>(y));
    float weight = Weight(d);
    auto weightgrad = WeightGrad(d);
    p.Velocity += g.Velocity * weight;
    p.Acceleration += g.Acceleration * weight;
    p.StrainRate[0][0] += g.Velocity.x * weightgrad.x;
    p.StrainRate[1][1] += g.Velocity.y * weightgrad.y;
    p.StrainRate[1][0] += (g.Velocity.x * weightgrad.y) + (g.Velocity.y * weightgrad.x);
    p.StrainRate[0][1] = p.StrainRate[1][0];
}
void G2P()
{
    for(int i = 0; i < ParticleCount;++i)
    {
        auto & particle = ParticleList[i];
        int GridX = std::floor(particle.Position.x);
        int GridY = std::floor(particle.Position.y);
        particle.Velocity = glm::vec2(0);
        particle.Acceleration = glm::vec2(0);
        particle.StrainRate = glm::mat2x2(0);
        for(int x = 0;x < 2;++x)
        {
            for(int y = 0;y < 2;++y)
            {
                if(GridX+x >=0 && GridX+x < GridSize)
                {
                    if(GridY+y >=0 && GridY+y < GridSize)
                    {
                        G2PNode(GridX+x,GridY+y,particle);
                    }
                }
            }
        }
    }
}
void Intergrate()
{
    for(int i = 0; i < ParticleCount;++i)
    {
        auto & particle = ParticleList[i];
//        std::cout << "i:p:"<<particle.Position.x<<" " << particle.Position.y<<std::endl;
//        std::cout << "i:v:"<<particle.Velocity.x<<" " << particle.Velocity.y<<std::endl;
        particle.Position += particle.Velocity * DeltaTime;
        particle.Velocity += particle.Acceleration * DeltaTime;
        particle.Strain += particle.StrainRate * DeltaTime;
        particle.Stress[0][0] = particle.Strain[0][0] * particle.E;
        particle.Stress[1][1] = particle.Strain[1][1] * particle.E;
        particle.Force = glm::vec2(0,-9.8 * particle.Mass);
    }
}
void Update()
{
    std::fill(SimGrid.begin(),SimGrid.end(),Grid());
    P2G();
    G2P();
    Intergrate();
}
void AddParticle(glm::vec2 pos){
    auto & p = ParticleList[ParticleCount];
    p.Position = pos;
    ParticleCount += 1;
}
void SaveParticles(GifWriter & g,int delay){

        for(int i = 0;i < ParticleCount;++i){
            auto p = ParticleList[i];
            for(int dx = 0;dx < 10;++dx){
                for(int dy = 0;dy < 10;++dy){
                    int x = static_cast<int>((p.Position.x) * Resolution) + dx;
                    int y = static_cast<int>((GridSize - p.Position.y) * Resolution) + dy;
                    int v = 4 * (x + (RenderSize * y));
                    if(x >= 0 && x < RenderSize){
                        if(y >= 0 && y < RenderSize){
                            frame[v] = 255;
                        }
                    }
                }
            }
        }
}
void PaintNumbers(int MaxTime,int t)
{
        float FillPC = (float)t / (float)MaxTime;
        static constexpr int Bits = 5;
        int fillint = static_cast<int>(FillPC * (float)(1<<(Bits-1)));
        std::bitset<Bits> pattern = std::bitset<Bits>(fillint);
        for(int b = 0;b < Bits;++b){
            for(int dx = 0;dx < 10;++dx){
                for(int dy = 0;dy < 10;++dy){
                    int x = static_cast<int>(b * Resolution) + dx;
                    int y = static_cast<int>(0) + dy;
                    int v = 4 * (x + (RenderSize * y));
                    if(x >= 0 && x < RenderSize){
                        if(y >= 0 && y < RenderSize){
                            frame[v+1] = 255;
                        }
                    }
                }
            }
        }
}
int main(int argc, char ** args)
{

	auto fileName = "out.gif";
	int delay = 10;
    for(int v = 0;v < 10;++v){
        AddParticle(glm::vec2(30.0+(v/2.0),40));
    }
    int MaxTime = 100; 
    GifWriter g;
	GifBegin(&g, fileName, RenderSize, RenderSize, delay);
    for(int t = 0;t < MaxTime;++t)
    {
        std::cout<<"T:"<<t<<std::endl;
        for(int i = 0; i< SubSteps;++i){
            Update();
        }
        std::fill(frame.begin(),frame.end(),0);
        SaveParticles(g,delay);
        PaintNumbers(MaxTime,t);
        GifWriteFrame(&g, frame.data(), GridSize * 10, GridSize * 10, delay);
    }
    std::fill(frame.begin(),frame.end(),255);
    GifWriteFrame(&g, frame.data(), GridSize * 10, GridSize * 10, delay);
	GifEnd(&g);
    std::cout<<"Finished"<<std::endl;
}
