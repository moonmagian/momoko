#include "pksystem.hpp"

momoko::pks::pksystem::pksystem(base::ideal_lattice &_latt)
    : latt(_latt), rng(std::random_device{}()) {}
