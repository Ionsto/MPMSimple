include <iostream>
#include "particle.h"
#include "grid.h"
#include "gif.h"
#include "SwapList.h"
#include <chrono>
#include <array>
#include <bitset>
#include <vector>
constexpr static int MaxParticles = 50000;
constexpr static int RealSize = 30;
constexpr static float GridDim = 0.5;
constexpr static int GridSize = static_cast<int>(static_cast<float>(RealSize)/GridDim);
constexpr static int MaxTime = 500;
//constexpr static float GridDim = static_cast<float>(RealSize) / static_cast<float>(GridSize);
constexpr static float DeltaTime = 1e-3;
static constexpr int SubSteps = 100;
static constexpr int Resolution = 20;
static constexpr int RenderSize = RealSize * Resolution;
//std::array<Particle,MaxParticles> ParticleList;
SwapList<Particle,MaxParticles> ParticleList;
std::array<Grid,GridSize * GridSize> SimGrid;
std::vector<uint8_t> frame(RenderSize * RenderSize * 4, 0);
Grid & GetGrid(int x, int y){
    return SimGrid[(x * GridSize) + y];
}
float WeightAxis(float x)
{
    return 1.0 - std::abs(x/GridDim);
}
float WeightGradAxis(float x)
{
    if(x==0)
    {
        return 0;
    }
    return std::copysign(-1.0,x) * GridDim;
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
    auto d = p.Position - glm::vec2(static_cast<float>(x) * GridDim,static_cast<float>(y) * GridDim);
    float weight = Weight(d);
    auto weightgrad = WeightGrad(d);
    g.Mass += p.Mass * weight;
    g.Velocity += p.Velocity * p.Mass * weight;
    g.Force += (p.Stress * weightgrad);
    g.Force += p.Force * weight;
}
void P2G()
{
    for(int i = 0; i < ParticleList.ParticleCount;++i)
    {
        auto & particle = ParticleList.Get(i);
        if(particle.Position.x < 0 || particle.Position.x > RealSize || particle.Position.y < 0 || particle.Position.y > RealSize)
        {
            ParticleList.Remove(i);
            i -= 1;
            continue;
        }
//        if(particle.Position.x < 0)
//        {
//            particle.Position.x = 0;
//        }
//        if(particle.Position.x >= GridSize)
//        {
//            particle.Position.x = GridSize - 0.1;
//        }
//        if(particle.Position.y < 0)
//        {
//            particle.Position.y = 0;
//        }
//        if(particle.Position.y >= GridSize)
//        {
//            particle.Position.y = GridSize - 0.1;
//        }
        int GridX = std::floor(particle.Position.x / GridDim);
        int GridY = std::floor(particle.Position.y / GridDim);
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
    float Friction = 0.1;
    for(int x = 0; x < GridSize;++x){
        for(int d = 0;d < 1;++d){
            //Floor 
            GetGrid(x,d).Velocity.y = 0;
            GetGrid(x,d).Velocity.x *= Friction;
            //Ceiling
            GetGrid(x,GridSize - (1 + d)).Velocity.y = 0;
            GetGrid(x,GridSize - (1 + d)).Velocity.x *= Friction;
            //Left wall
            GetGrid(d,x).Velocity.x = 0; 
            GetGrid(d,x).Velocity.y *= Friction;
            //Right wall
            GetGrid(GridSize - (1 + d),x).Velocity.x = 0;
            GetGrid(GridSize - (1 + d),x).Velocity.y *= Friction;
        }
    }
}
void G2PNode(int x,int y,Particle & p){
    Grid g = GetGrid(x,y);
    auto d = p.Position - glm::vec2(static_cast<float>(x * GridDim),static_cast<float>(y * GridDim));
    float weight = Weight(d);
    auto weightgrad = WeightGrad(d);
    p.Velocity += g.Velocity * weight;
    p.Acceleration += g.Acceleration * weight;
//    p.Acceleration.x += -1e3 * g.Mass * weightgrad.x;
//    p.Acceleration.y += -1e3 * g.Mass * weightgrad.y;
    p.StrainRate[0][0] += g.Velocity.x * weightgrad.x;
    p.StrainRate[1][1] += g.Velocity.y * weightgrad.y;
    p.StrainRate[1][0] += 0.5 * (g.Velocity.x * weightgrad.y) + (g.Velocity.y * weightgrad.x);
    p.StrainRate[0][1] = p.StrainRate[1][0];
}
void G2P()
{
    for(int i = 0; i < ParticleList.ParticleCount;++i)
    {
        auto & particle = ParticleList.Get(i);
        int GridX = std::floor(particle.Position.x / GridDim);
        int GridY = std::floor(particle.Position.y / GridDim);
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
    for(int i = 0; i < ParticleList.ParticleCount;++i)
    {
        auto & particle = ParticleList.Get(i);
//        std::cout << "i:p:"<<particle.Position.x<<" " << particle.Position.y<<std::endl;
//        std::cout << "i:v:"<<particle.Velocity.x<<" " << particle.Velocity.y<<std::endl;
        particle.Position += particle.Velocity * DeltaTime;
        particle.Velocity += particle.Acceleration * DeltaTime;
        particle.Strain += particle.StrainRate * DeltaTime;
        particle.Force = glm::vec2(0,-9.8 * particle.Mass);
//Calculate stress
        if(particle.Type == 0){
            float prefix = (particle.YoungsModulus) / (1 - (particle.PoissonsRatio * particle.PoissonsRatio));
            particle.Stress[0][0] = prefix * (particle.Strain[0][0] + (particle.PoissonsRatio * particle.Strain[1][1]));
            particle.Stress[1][1] = prefix * (particle.Strain[1][1] + (particle.PoissonsRatio * particle.Strain[0][0]));
            float meanshearstrain = (particle.Strain[0][1] + particle.Strain[1][0]) / 2.0;
            float deltaxy = ((0.5 * particle.YoungsModulus) / (1.0 + particle.PoissonsRatio)) * meanshearstrain;
            particle.Stress[0][1] = deltaxy;
            particle.Stress[1][0] = deltaxy;
        }
        else{
            float pressure = particle.YoungsModulus * 0.5 * (particle.Strain[0][0] + particle.Strain[1][1]);
            particle.Stress[0][0] = pressure;
            particle.Stress[1][1] = pressure;
            float meanshearstrainrate = (particle.StrainRate[0][1] + particle.StrainRate[1][0]) / 2.0;
            float shear = particle.Viscosity * meanshearstrainrate;
            particle.Stress[0][1] = shear;
            particle.Stress[1][0] = shear;
        }

//Check in bounds
//        if(particle.Position.x < 0)
//        {
//            particle.Position.x = 0;
//        }
//        if(particle.Position.x >= GridSize)
//        {
//            particle.Position.x = GridSize - 0.1;
//        }
//        if(particle.Position.y < 0)
//        {
//            particle.Position.y = 0;
//        }
//        if(particle.Position.y >= GridSize)
//        {
//            particle.Position.y = GridSize - 0.1;
//        }
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
    auto p = Particle();
    p.Position = pos;
    p.Type = 1; 
    ParticleList.Add(p);
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
void SaveParticles(GifWriter & g,int delay){

        for(int i = 0;i < ParticleList.ParticleCount;++i){
            auto p = ParticleList.Get(i);
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
int main(int argc, char ** args)
{

	auto fileName = "out.gif";
	int delay = 10;
    int count = 10;
    float density = 1.5;
//    for(int v = 0;v < count;++v){
//        for(int u = 0;u < count;++u){
//            AddParticle(glm::vec2(10.0+(u/density),30+(v/density)));
//        }
//    }
//    count = 5;
//    density = 1;
//    for(int v = 1;v < (GridSize - 1)* density;++v){
//        for(int u = 0;u < 5*density;++u){
//            AddParticle(glm::vec2((v/density),(u/density)));
//        }
//    }
    GifWriter g;
	GifBegin(&g, fileName, RenderSize, RenderSize, delay);
    using namespace std::chrono;
    auto start = high_resolution_clock::now();
    for(int t = 0;t < MaxTime;++t)
    {
        for(int p = 0;p < 80;++p)
        {
            auto pa = Particle();
            float dx = 2*((rand() % 1000) / 1000.0);
            float dy = 2*((rand() % 1000) / 1000.0);
            pa.Position = glm::vec2(10.0+dx,20.0+dy); 
            pa.Velocity.x = 0;
            pa.Type = 1;
            pa.Colour.r = rand()%255;
            pa.Colour.g = rand()%255;
            pa.Colour.b = rand()%255;
            ParticleList.Add(pa);
        }
        std::cout<<"T:"<<t<<std::endl;
        for(int i = 0; i < SubSteps;++i){
            Update();
        }
        std::fill(frame.begin(),frame.end(),0);
        SaveParticles(g,delay);
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
