{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}

import Prelude hiding (sum, concat, mapM, sequence)
import Data.Foldable
import Data.Traversable
import Data.Complex
import Control.Applicative
import Control.Lens
import Linear
import Linear.V (V)
import qualified Linear.V as LV
import qualified Data.Vector as V
import GHC.TypeLits (Nat)
import System.Random.MWC hiding (uniform)
import Data.Random hiding (uniform)

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

uniform :: RealFloat a => Models a
uniform = l1Normalize $ pure 1

noisify :: (Traversable f, Num a, Distribution Normal a)
        => a -> f (Point x a) -> IO (f (Point x a))
noisify sigma pts =
    createSystemRandom >>= runRVar (mapM addNoise pts)
  where
    addNoise p = (\a->p { pY=pY p + a}) <$> normal 0 sigma
  
main = do
    let f = uniform
    let f = l1Normalize $ v [2,2,3,2, 1,8]
    pts' <- noisify 0.1 pts
    writeFile "points.txt" $ unlines
      $ map (\(Point (V1 x) y s)->show x++"\t"++show y++"\t"++show s) pts'
    print -- $ takeWhile (\f->test pts' models f > 0.1)
          $ maxEnt 1e-2 pts' models f

testGrad = do
    let f :: RealFloat a => Models a
        f = l1Normalize $ pure 1
    let gradError n = 
          let pts = samplePts n
              fin = finiteGrad 1e-14 (chiSquared pts models) (f :: Models Double)
              hand = gradChiSquared pts models f
              ad = grad (chiSquared (fmap (fmap realToFrac) pts) models) f
              put msg a = putStrLn $ msg++"\t"++show a
              relErr a b = norm (a ^-^ b) / norm b
          in do put "fin-hand" $ relErr fin hand
                put "fin-ad"   $ relErr fin ad
                put "ad-hand"  $ relErr ad hand
    Data.Foldable.mapM_ gradError [2, 10, 20, 50, 100, 200, 400, 800, 2000, 4000, 8000]

testHessian = do
    let f :: RealFloat a => Models a
        f = l1Normalize $ pure 1
    let gradError n = 
          let pts = samplePts n
              hand = hessianChiSquared pts models f
              ad = hessian (chiSquared (fmap (fmap realToFrac) pts) models) f
              put msg a = putStrLn $ msg++"\t"++show a
              frobenius = sum . fmap (sum . (^2))
              relErr a b = frobenius (a !-! b) / frobenius b
          in put "ad-hand"  $ relErr ad hand
    Data.Foldable.mapM_ gradError [2, 10, 20, 50, 100, 200, 400, 800, 2000, 4000, 8000]
    
