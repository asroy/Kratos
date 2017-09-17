#pragma once
#include <stdio.h>
#include <string>

namespace OversetAssembly
{
struct DonorInfo
{
    std::size_t mNumDonorNode;
    std::vector<std::size_t> mDonorNodesEquationId;
    double mBarycentricCoordinate[3];
};

}