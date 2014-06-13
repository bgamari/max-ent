{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE FlexibleInstances #-}

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
import FiniteDiff
import Data.Number.BigFloat
import Numeric.AD

instance Linear.Epsilon (BigFloat Prec50) where nearZero = (<1e-14)

expModel :: RealFloat a => a -> V1 a -> a
expModel tau (V1 t) = exp (-t/tau)

v :: LV.Dim n => [a] -> V (n::Nat) a
v = maybe (error "Invalid length vector given") id . LV.fromVector . V.fromList

type Models = V 6

models :: Models (Model V1)
models = v [ Model $ expModel 100, Model $ expModel 300
           , Model $ expModel 500, Model $ expModel 700
           , Model $ expModel 900, Model $ expModel 1100
           ]
           
weights :: RealFloat a => Models a
weights = l1Normalize $ v [1,2,10,2, 1,9]

samplePts n =
    map (\x->Point (V1 x) (mixtureModel models weights (V1 x)) 1)
    $ linSpace 0 2000 n

pts = samplePts 2000

linSpace :: RealFrac a => a -> a -> Int -> [a]
linSpace a b n = [a + (b-a)/realToFrac n*realToFrac i | i <- [0..n]]

main = do
    let f :: RealFloat a => Models a
        f = l1Normalize $ pure 1

    --putStrLn $ unlines $ map (\(Point (V1 x) y s)->show x++"\t"++show y) pts

    print $ takeWhile (\f->test pts models f > 0.1)
          $ maxEnt 1 pts models f

testGrad = do
    let f :: RealFloat a => Models a
        f = l1Normalize $ pure 1
    let gradError n = 
          let pts = samplePts n
              fin = finiteGrad 1e-14 (chiSquared pts models) (f :: Models Double)
              hand = gradChiSquared pts models f
              ad = grad (chiSquared (fmap (fmap realToFrac) pts) models) f
              put msg a = putStrLn $ msg++"\t"++show a
          in do put "fin-hand" $ norm (fin ^-^ hand) / norm fin
                put "fin-ad"   $ norm (fin ^-^ ad) / norm fin
                put "ad-hand"  $ norm (hand ^-^ ad) / norm fin
    Data.Foldable.mapM_ gradError [2, 10, 20, 50, 100, 200, 400, 800, 2000, 4000, 8000]
