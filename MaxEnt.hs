{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ImpredicativeTypes #-}

module MaxEnt where

import Prelude hiding (sum)
import Data.Foldable as F
import Data.Traversable as T
import Data.Monoid
import Control.Applicative
import Linear

import BunseGurstner

data Point x a = Point { pX :: !(x a)
                       , pY, pSigma :: !a
                       }

instance Functor x => Functor (Point x) where
  fmap f (Point x y s) = Point (fmap f x) (f y) (f s)

-- | @y = model x@
data Model x = Model {runModel :: forall a. RealFloat a => x a -> a}

onValue :: RealFloat a => x a -> Model x -> a
onValue = flip runModel

mixtureModel :: (Foldable f, Applicative f, RealFloat a)
             => f (Model x)   -- ^ models to mix
             -> f a           -- ^ mixing weights
             -> x a           -- ^ abscissa
             -> a             -- ^ ordinate
mixtureModel models fs x = sum $ (\f m->f * runModel m x) <$> fs <*> models

sumM :: (Num a, Foldable f, Applicative v, Additive w, Additive v)
     => f (v (w a)) -> v (w a)
sumM = F.foldl' (!+!) (pure zero)


entropy :: (Foldable f, RealFloat a) => f a -> a
entropy = getSum . F.foldMap (\p->Sum $ p * log p)

gradEntropy :: (Functor f, RealFloat a) => f a -> f a
gradEntropy = fmap (\x->1 + log x)
            
hessianEntropy :: (Traversable f, RealFloat a) => f a -> f (f a)
hessianEntropy = kronecker . fmap recip               


chiSquared :: (Foldable f, Applicative f, RealFloat a)
           => [Point x a]     -- ^ points
           -> f (Model x)     -- ^ models
           -> f a             -- ^ weights
           -> a
chiSquared pts models f = getSum $ foldMap g pts
  where
    g (Point x y s) = Sum $ (mixtureModel models f x - y)^2 / s^2

gradChiSquared :: (Traversable f, Additive f, Applicative f, RealFloat a)
               => [Point x a] -> f (Model x) -> f a -> f a
gradChiSquared pts models f = (-2) *^ sumV (fmap g pts)
  where
    g (Point x y s) = (mixtureModel models f x + y) / s *^ fmap (onValue x) models

hessianChiSquared :: (RealFloat a, Additive f, Applicative f)
                  => [Point x a] -> f (Model x) -> f a -> f (f a)
hessianChiSquared pts models f = (-2) *!! sumM (fmap g pts)
  where
    g (Point x y s) =
      let m = fmap (onValue x) models
      in (m `outer` m) !!* recip s


-- | Maximum entropy fitting of distribution
maxEnt :: (Metric f, Traversable f, Applicative f, RealFloat a)
       => [Point x a]           -- ^ observations
       -> f (Model x)           -- ^ models
       -> f a                   -- ^ initial weights
       -> [f a]                 -- ^ iterates
maxEnt pts models = go
  where
    go f = undefined -- f ^+^ sumV (dot <$> ZipList basis <*> x)
      where
        basis = toList $ subspace pts f models
        x = ZipList undefined

-- | Component-wise multiplication
cmult :: (Num a, Applicative f) => f a -> f a -> f a
cmult a b = (*) <$> a <*> b

cmultOp :: (Num a, Applicative f) => f a -> f (f a) -> f (f a)
cmultOp v a = (*^) <$> v <*> a

subspace :: (Applicative f, Traversable f, Metric f, RealFloat a)
         => [Point x a] -> f a -> f (Model x) -> V3 (f a)
subspace pts f models = V3 e1 e2 e3
  where
    e1 = f `cmult` gs
    e2 = f `cmult` gc
    e3 =   d gs *^ ((f `cmultOp` hc) !* (f `cmult` gs))
       ^-^ d gc *^ ((f `cmultOp` hc) !* (f `cmult` gc))
    hc = hessianChiSquared pts models f
    gs = gradEntropy f
    gc = gradChiSquared pts models f
    d grad = recip $ sqrt $ F.sum $ f `cmult` fmap (^2) grad

-- | Eigenvalues of 3x3 symmetric matrix 
eigen3 :: (Epsilon a, RealFloat a) => V3 (V3 a) -> [a]
eigen3 a = cubicTrig (-3) (det33 a)

-- | Project the given operator into a new basis
projectInto :: (Foldable g, Metric g,
                Foldable f, Metric f,
                Num a)
            => g (g a) -> f (g a) -> f (f a)
projectInto op basis = fmap (\em->fmap (f em) basis) basis
  where
    f em ev = em `dot` (op !* ev)
