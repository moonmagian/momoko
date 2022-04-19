#ifndef DESCRETE_GAUSSIAN_CALCULATOR_HPP
#define DESCRETE_GAUSSIAN_CALCULATOR_HPP

namespace momoko::gaussian {}
class descrete_gaussian_calculator {
private:
  double delta;
  double tail;

public:
  descrete_gaussian_calculator(double delta, double tail);
};

#endif // DESCRETE_GAUSSIAN_CALCULATOR_HPP
