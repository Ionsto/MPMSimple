#if(
#{
#    return 0.125 * std::pow(3 + ((2*x)/GridDim),2);
#}
#if(x > (GridDim)/2.0 && x < (3.0f * GridDim)/2.0)
#{
#    return 0.125 * std::pow(3 - ((2*x)/GridDim),2);
#}
#if(x >= -(GridDim)/2.0 && x <= (GridDim)/2.0)
#{
#    return 0.75 - std::pow((x/GridDim),2);
#}
#return 0;
import matplotlib.pyplot as plt
import numpy as np
def shape(x):
    GridDim = 1
    if (x < (-(GridDim)/2.0)) and (x > -(3.0 * GridDim)/2.0):
        return 0.125 * np.power(3 + ((2*x)/GridDim),2)
    if (x > ((GridDim)/2.0)) and (x < (3.0 * GridDim)/2.0):
        return 0.125 * np.power(3 - ((2*x)/GridDim),2)
    if (x >= -(GridDim)/2.0) and (x <= (GridDim)/2.0):
        return 0.75 - np.power((x/GridDim),2)
    return 0
