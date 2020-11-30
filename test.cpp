#include <iostream>
#include <math.h>
float GridDim = 1;
float WeightAxis(float x)
{
    if(x < -(GridDim)/2.0 && x > -(3.0f * GridDim)/2.0)
    {
        return 0.125 * std::pow(3 + ((2*x)/GridDim),2);
    }
    if(x > (GridDim)/2.0 && x < (3.0f * GridDim)/2.0)
    {
        return 0.125 * std::pow(3 - ((2*x)/GridDim),2);
    }
    if(x >= -(GridDim)/2.0 && x <= (GridDim)/2.0)
    {
        return 0.75 - std::pow((x/GridDim),2);
    }
    return 0;
}
int main(int argc,char **args){
    while(true){
        float i;
        std::cin>>i;
        std::cout<<WeightAxis(i)<<std::endl;
    }
}
