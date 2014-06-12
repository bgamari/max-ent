import Prelude hiding (sum)
import Data.Foldable
import Data.Traversable
import Data.Complex
import Control.Applicative
import Linear

import MaxEnt

expModel :: RealFloat a => a -> V1 a -> a
expModel tau (V1 t) = exp (-t/tau)

models = V4 (Model $ expModel 100) (Model $ expModel 200)
            (Model $ expModel 500) (Model $ expModel 700)

weights :: V4 Double
weights = normalize $ V4 1 1 1 1
pts = map (\x->Point (V1 x) (mixtureModel models weights (V1 x)) 1) [0..2000]

main = do
    let f = normalize $ V4 2 3 2 3
    let hc = hessianChiSquared pts models f
        basis = subspace pts f models
    print $ gradChiSquared pts models f ^-^ finiteGrad 1e-7 (chiSquared pts models) f
    print $ gradChiSquared pts models f

    print $ maxEnt 1e3 pts models f


finiteGrad :: (Traversable f, Applicative f, Additive f, RealFrac a, Epsilon a)
           => a -> (f a -> a) -> f a -> f a
finiteGrad tol f x = fmap (finiteDiff tol f x) $ kronecker (pure 1)

finiteDiff :: (Traversable f, Applicative f, Additive f, RealFrac a, Epsilon a)
           => a -> (f a -> a) -> f a -> f a -> a
finiteDiff tol f x dir = go 1
  where
    alpha = 2
    diff h = (f (x ^+^ h*^dir) - f x) / h
    go h
      | nearZero h = 1/0
      | err < tol  = d1
      | otherwise  = go (h/alpha)
      where
        d0 = diff h
        d1 = diff (h*alpha)
        err = abs (d0 - d1) / d1
