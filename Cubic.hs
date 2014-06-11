module Cubic (cubicTrig) where

import Linear

-- | Trignometric solutions to a cubic equation `x^3 + p*x + q == 0`
cubicTrig :: (Epsilon a, RealFloat a) => a -> a -> [a]
cubicTrig p q
  -- x^3 == 0
  | nearZero q && nearZero p = [0]
  -- x^3 + p*x == 0
  | nearZero (q^2) && p > 0 = []
  | nearZero (q^2) = [sqrt (-p)]
  -- x^3 + q == 0
  | nearZero (p^2) && q > 0 = []
  | nearZero (p^2) = [(-q)**(1/3)]
  -- x^3 + p*x + q == 0
  | 4 * p^3 + 27*q^2 <= 0 =
    let t k = 2 * sqrt (-p/3) * cos (theta/3 - 2*pi*k/3)
        theta = acos (3/2 * q/p * sqrt (-3 / p))
    in map t [0,1,2]
  | p > 0 = 
    let theta = asinh (3/2 * q/p * sqrt (3 / p))
    in [-2 * sqrt (p/3) * sinh (theta/3)]
  | otherwise =
    let theta = acosh (-3/2 * abs q / p * sqrt (-3/p))
    in [-2 * signum q * sqrt (-p/3) * cosh (theta / 3)]

