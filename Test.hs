{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}

import Prelude hiding (sum, concat)
import Data.Foldable
import Data.Traversable
import Data.Complex
import Control.Applicative
import Linear
import Linear.V (V)
import qualified Linear.V as LV
import qualified Data.Vector as V
import GHC.TypeLits (Nat)

import MaxEnt

expModel :: RealFloat a => a -> V1 a -> a
expModel tau (V1 t) = exp (-t/tau)

v :: LV.Dim n => [a] -> V (n::Nat) a
v = maybe (error "Invalid length vector given") id . LV.fromVector . V.fromList

models :: V 8 (Model V1)
models = v [ Model $ expModel 100, Model $ expModel 200
           , Model $ expModel 300, Model $ expModel 400
           , Model $ expModel 500, Model $ expModel 600
           , Model $ expModel 700, Model $ expModel 800
           ]

weights :: V 8 Double
weights = l1Normalize $ v [1,2,1,10, 1,2,8,1]

pts = map (\x->Point (V1 x) (mixtureModel models weights (V1 x)) 1) [0..2000]

main = do
    let f :: V 8 Double
        f = l1Normalize $ v $ concat $ replicate 8 [1]
    let hc = hessianChiSquared pts models f
        basis = subspace pts f models
    --print $ gradChiSquared pts models f ^-^ finiteGrad 1e-7 (chiSquared pts models) f

    putStrLn $ unlines $ map (\(Point (V1 x) y s)->show x++"\t"++show y) pts
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
