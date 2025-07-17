#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>




namespace CPS {
namespace EMT {
namespace Ph1 {


class PieceWiseNonLinearCharacteristic {
public:
    // A point is given as (current, flux linkage).
    // Points must be provided in order of increasing current.
    struct Point {
        double current; // x-axis (per unit)
        double flux;    // y-axis (per unit)
    };

    // Constructor accepts a vector of boundary points and a saturated inductance value.
    // The points define the boundaries for the piecewise linear region.
    // If the flux is higher than the last boundary point, the saturated region is used.
    PieceWiseNonLinearCharacteristic(const std::vector<Point>& pts, double satL);

    // Compute the inductance for a given flux linkage.
    // For flux values within an interval: L = (flux₂ - flux₁)/(current₂ - current₁).
    // For flux values above the last point, the saturated inductance is returned.
    double getInductance(double flux) const;

    // Compute the current for a given flux linkage.
    // For flux values within an interval, the current is computed by linear interpolation.
    // For flux values above the last boundary, the current is:
    // I = I_last + (flux - flux_last) / mSaturatedInductance.
    // If flux is negative, the computed current is returned as negative.
    double getCurrent(double flux) const;

private:
    std::vector<Point> mPoints;
    double mSaturatedInductance;
};


} // namespace Ph1
} // namespace EMT
} // namespace CPS