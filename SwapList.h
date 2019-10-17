#pragma once
#include <iostream>
#include <array>
template<class t, int max>
class SwapList {
public:
	static constexpr int MaxParticlesPerElement = max;
	std::array<t, MaxParticlesPerElement> ParticleList;
	int ParticleCount = 0;
	t& Get(int i)
	{
		return ParticleList[i];
	}
	std::array<t, MaxParticlesPerElement> RawData()
	{
		return ParticleList;
	}
	bool Add(t add)
	{
		if (ParticleCount != MaxParticlesPerElement)
		{
			ParticleList[ParticleCount++] = add;
			return true;
		}
		std::cout << "Failed to add particle\n";
		return false;
	}
	void Remove(int position)
	{
		if (position == ParticleCount - 1)
		{
			ParticleCount--;
			return;
		}
		if (position < ParticleCount)
		{
			ParticleList[position] = ParticleList[ParticleCount-- - 1];
		}
		else {
			throw;
			//wtf
		}
	}
};
