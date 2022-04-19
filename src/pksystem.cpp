#include "pksystem.hpp"

momoko::pks::pksystem::pksystem(base::ideal_lattice &_latt,
                                gaussian::gaussian_dist_sampler &_sampler)
    : latt(_latt), sampler(_sampler), rng(std::random_device{}()) {}
