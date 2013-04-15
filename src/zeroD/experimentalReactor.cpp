#include "cantera/zeroD/ExperimentalReactor.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/thermo/SurfPhase.h"

#include <cfloat>

using namespace std;

namespace Cantera {

ExperimentalReactor::ExperimentalReactor()
{
    cout << "Creating Experimental Reactor" << endl;
    m_mass = -1;
}

void ExperimentalReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total internal energy,
    // [1] the total volume, and [2...K+2] the mass of each species.
    m_vol = y[1];

    // Set the mass fractions
    if (m_mass < 0) {
        m_mass = accumulate(y+2, y+2+m_nsp, 0.0);
    }

    vector_fp yNow(y+2, y+2+m_nsp);
    for (size_t k = 0; k < m_nsp; k++) {
        yNow[k] /= m_mass;
    }
    m_thermo->setMassFractions_NoNorm(&yNow[0]);

    if (m_energy) {
        // Use Newton's method to determine the mixture temperature. Tight
        // tolerances are required both for Jacobian evaluation and for
        // sensitivity analysis to work correctly.

        doublereal U = y[0];
        doublereal T = temperature();
        double dT = 100;

        int i = 0;
        while (abs(dT / T) > 10 * DBL_EPSILON) {
            m_thermo->setState_TR(T, m_mass / m_vol);
            double dUdT = m_thermo->cv_mass() * m_mass;
            dT = (m_thermo->intEnergy_mass() * m_mass - U) / dUdT;
            T -= dT;
            i++;
            if (i > 100) {
                throw CanteraError("Reactor::updateState", "no convergence");
            }
        }
    } else {
        m_thermo->setDensity(m_mass/m_vol);
    }

    size_t loc = m_nsp + 2;
    SurfPhase* surf;
    for (size_t m = 0; m < m_nwalls; m++) {
        surf = m_wall[m]->surface(m_lr[m]);
        if (surf) {
            m_wall[m]->setCoverages(m_lr[m], y+loc);
            loc += surf->nSpecies();
        }
    }

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}

}
