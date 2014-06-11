import Prelude hiding (sum)                
import Data.Foldable
import Data.Complex
import Control.Applicative
import Linear

import Test.QuickCheck.Property
import Test.QuickCheck.Modifiers

import MaxEnt

expModel :: RealFloat a => a -> V1 a -> a
expModel tau (V1 t) = exp (-t/tau)         

models = V4 (Model $ expModel 100) (Model $ expModel 200)
            (Model $ expModel 500) (Model $ expModel 700)
weights :: V4 Double
weights = normalize $ V4 1 1 1 1
pts = map (\x->Point (V1 x) (mixtureModel models weights (V1 x)) 1) [0..2000]

main = do
    let f = weights
    let hc = hessianChiSquared pts models f
        basis = subspace pts f models
        m = hc `projectInto` basis
        g = fmap (\em->fmap (em `dot`) basis) basis
    print $ eigen3 g
    print $ eigen3 m
    
cubicSolves :: Double -> Double -> Property
cubicSolves p q
  | nearZero p && nearZero q = property rejected
  | otherwise =
    conjoin $ map solves $ cubicTrig p q
  where
     cubic x = x^3 + p*x + q
     solves :: Double -> Property
     --solves x | isNaN x = property rejected
     solves x = counterexample msg $ nearZero $ (cubic x)^2
       where msg = show x++" is not a root (f(x) = "++show (cubic x)++")"

runTests = [property cubicSolves]
