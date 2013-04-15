#include "cantera/zeroD/Reactor.h"

namespace Cantera {

class ExperimentalReactor : public Reactor
{
public:
    ExperimentalReactor();
    virtual int type() const { return ExperimentalReactorType; }

    virtual void updateState(doublereal* y);

private:
    double m_mass;
};

}
