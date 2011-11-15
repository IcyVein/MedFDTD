#include "computeFraction.h"

double computeFraction(const double *PartMass, const double requiredMass)
{
    double massThreshold = 0.01, currentMass = 0, deltaMass = 10;
    double fraction = 0.5, fractionTemp = fraction, fractionStep = fraction/2;

    while (deltaMass > massThreshold)
    {
        fraction = fractionTemp;
        currentMass = PartMass[3]*fraction*fraction*fraction+PartMass[2]*fraction*fraction+PartMass[1]*fraction+PartMass[0];
        if (currentMass > requiredMass)
            fractionTemp = fraction-fractionStep;
        else if (currentMass < requiredMass)
            fractionTemp = fraction+fractionStep;
        else
            break;
        fractionStep = fractionStep/2;
        if (currentMass > requiredMass)
            deltaMass = currentMass-requiredMass;
        else if (currentMass < requiredMass)
            deltaMass = requiredMass-currentMass;
        else
            deltaMass = 0;
    }
    return fraction;
}

