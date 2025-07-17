#include "dpsim-models/EMT/EMT_Ph1_PieceWiseNonLinearCharacteristic.h"

// Constructor: Validates that the points are nonempty and that current increases.
CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::PieceWiseNonLinearCharacteristic(const std::vector<Point>& pts, double satL)
    : mPoints(pts), mSaturatedInductance(satL)
{
    if(mPoints.empty())
        throw std::invalid_argument("At least one point must be provided.");
    // Check that currents increase.
    for (size_t i = 1; i < mPoints.size(); ++i) {
        if (mPoints[i].current <= mPoints[i - 1].current)
            throw std::invalid_argument("Points must have strictly increasing currents.");
    }
}

double CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::getInductance(double flux) const {
    // Work with the absolute flux for determination, then mirror later if needed.
    bool isNegative = (flux < 0);
    double absFlux = std::fabs(flux);

    // Search through the defined intervals.
    size_t numIntervals = mPoints.size() - 1; // There are n-1 intervals for n points.
    for (size_t i = 0; i < numIntervals; ++i) {
        if (absFlux <= mPoints[i + 1].flux) {
            double deltaFlux = mPoints[i + 1].flux - mPoints[i].flux;
            double deltaCurrent = mPoints[i + 1].current - mPoints[i].current;
            if (deltaCurrent == 0)
                throw std::runtime_error("Zero current difference in interval - cannot compute inductance.");
            return deltaFlux / deltaCurrent;
        }
    }
    // Flux is higher than all defined boundaries - use saturated inductance.
    return mSaturatedInductance;
}

double CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::getCurrent(double flux) const {
    bool isNegative = (flux < 0);
    double absFlux = std::fabs(flux);
    double current = 0.0;
    
    size_t numIntervals = mPoints.size() - 1;
    for (size_t i = 0; i < numIntervals; ++i) {
        if (absFlux <= mPoints[i + 1].flux) {
            // Perform linear interpolation in this interval.
            double deltaFlux = mPoints[i + 1].flux - mPoints[i].flux;
            if (deltaFlux == 0)
                throw std::runtime_error("Zero flux difference in interval - cannot compute interpolation.");
            double t = (absFlux - mPoints[i].flux) / deltaFlux;
            current = mPoints[i].current + t * (mPoints[i + 1].current - mPoints[i].current);
            return isNegative ? -current : current;
        }
    }
    // AbsFlux is above the last boundary: we are in saturation.
    // Compute current by extending the last segment using the saturated inductance.
    double lastCurrent = mPoints.back().current;
    double lastFlux = mPoints.back().flux;
    current = lastCurrent + (absFlux - lastFlux) / mSaturatedInductance;
    return isNegative ? -current : current;
}