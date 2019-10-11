#include <iostream>
#include "particle.h"
#include "grid.h"
#include <array>
constexpr static int MaxParticles = 1000;
constexpr static int GridSize = 100;
constexpr static float DeltaTime = 1e-2;
int ParticleCount = 0;
std::array<Particle,MaxParticles> ParticleList;
std::array<Grid,GridSize * GridSize> SimGrid;
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
void SaveParticles(){

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
        }
    }
    for(int x = 0; x < GridSize;++x){
        GetGrid(x,1).Velocity.y = 0;
    }
}
void G2PNode(int x,int y,Particle p){
    Grid & g = GetGrid(x,y);
    auto d = p.Position - glm::vec2(static_cast<float>(x),static_cast<float>(y));
    float weight = Weight(d);
    auto weightgrad = WeightGrad(d);
    p.Velocity += g.Velocity * weight;
    p.StrainRate[0][0] += g.Velocity.x * weightgrad.x;
    p.StrainRate[1][1] += g.Velocity.y * weightgrad.y;
    p.StrainRate[1][0] += (g.Velocity.x * weightgrad.y) + (g.Velocity.y * weightgrad.x);
    p.StrainRate[0][1] = p.StrainRate[1][0];
}
void G2P()
{
    for(int i = 0; i < ParticleCount;++i)
    {
        auto particle = ParticleList[i];
        int GridX = std::floor(particle.Position.x);
        int GridY = std::floor(particle.Position.y);
        particle.Velocity = glm::vec2();
        particle.StrainRate = glm::mat2x2();
        for(int x = 0;x < 2;++x)
        {
            for(int y = 0;y < 2;++y)
            {
                G2PNode(GridX+x,GridY+y,particle);
            }
        }
    }
}
void Intergrate()
{
    for(int i = 0; i < ParticleCount;++i)
    {
        auto particle = ParticleList[i];
        particle.Position += particle.Velocity * DeltaTime;
        particle.Velocity += particle.Acceleration * DeltaTime;
        particle.Strain += particle.StrainRate * DeltaTime;
        particle.Stress[0][0] = particle.Strain[0][0] * particle.E;
        particle.Stress[1][1] = particle.Strain[1][1] * particle.E;
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
int main(int argc, char ** args)
{
    for(int v = 0;v < 10;++v){
        AddParticle(glm::vec2(30+(v/2.0),40));
    }
    int MaxTime = 100; 
    for(int t = 0;t < MaxTime;++t)
    {
        std::cout<<"T:"<<t<<std::endl;
        Update();
        SaveParticles();
    }
    std::cout<<"Finished"<<std::endl;
}
