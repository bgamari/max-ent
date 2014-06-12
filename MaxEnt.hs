{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ImpredicativeTypes #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}

module MaxEnt where

import Prelude hiding (sum, any)
import Data.Foldable as F
import Data.Traversable as T
import Data.Distributive
import Data.Monoid
import Control.Applicative
import Control.Lens
import Linear

import qualified Data.Packed.Matrix as HM
import qualified Data.Packed.Vector as HV
import qualified Numeric.LinearAlgebra.LAPACK as LA
import Foreign.Storable

import Debug.Trace
import Unsafe.Coerce

data Point x a = Point { pX :: !(x a)
                       , pY, pSigma :: !a
                       }

instance Functor x => Functor (Point x) where
  fmap f (Point x y s) = Point (fmap f x) (f y) (f s)

-- | @y = model x@
data Model x = Model {runModel :: forall a. RealFloat a => x a -> a}

onValue :: RealFloat a => x a -> Model x -> a
onValue = flip runModel
{-# INLINE onValue #-}

mixtureModel :: (Foldable f, Applicative f, RealFloat a)
             => f (Model x)   -- ^ models to mix
             -> f a           -- ^ mixing weights
             -> x a           -- ^ abscissa
             -> a             -- ^ ordinate
mixtureModel models fs x = sum $ (\f m->f * runModel m x) <$> fs <*> models
{-# INLINE mixtureModel #-}

sumM :: (Num a, Foldable f, Applicative v, Additive w, Additive v)
     => f (v (w a)) -> v (w a)
sumM = F.foldl' (!+!) (pure zero)
{-# INLINE sumM #-}

entropy :: (Foldable f, RealFloat a) => f a -> a
entropy = getSum . F.foldMap (\p->Sum $ p * log p)
{-# INLINE entropy #-}

gradEntropy :: (Functor f, RealFloat a) => f a -> f a
gradEntropy = fmap (\p->1 + log p)
{-# INLINE gradEntropy #-}

hessianEntropy :: (Traversable f, RealFloat a) => f a -> f (f a)
hessianEntropy = kronecker . fmap recip
{-# INLINE hessianEntropy #-}


chiSquared :: (Foldable f, Applicative f, RealFloat a)
           => [Point x a]     -- ^ points
           -> f (Model x)     -- ^ models
           -> f a             -- ^ weights
           -> a
chiSquared pts models f = getSum $ foldMap g pts
  where
    g (Point x y s) = Sum $ (mixtureModel models f x - y)^2 / s^2
{-# INLINE chiSquared #-}

gradChiSquared :: (Traversable f, Additive f, Applicative f, RealFloat a)
               => [Point x a] -> f (Model x) -> f a -> f a
gradChiSquared pts models f = 2 *^ sumV (fmap g pts)
  where
    g (Point x y s) = (mixtureModel models f x - y) / s *^ fmap (onValue x) models
{-# INLINE gradChiSquared #-}

hessianChiSquared :: (RealFloat a, Additive f, Applicative f)
                  => [Point x a] -> f (Model x) -> f a -> f (f a)
hessianChiSquared pts models f = 2 *!! sumM (fmap g pts)
  where
    g (Point x y s) =
      let m = fmap (onValue x) models
      in (m `outer` m) !!* recip s
{-# INLINE hessianChiSquared #-}


showMatrix :: (Functor f, Foldable f, Foldable g, Show a) => f (g a) -> String
showMatrix a = unlines $ toList $ fmap (foldMap (pad 70 . show)) a
  where pad n s = take n $ s ++ repeat ' '

-- | Maximum entropy fitting of distribution
maxEnt :: forall f x a.
          (Distributive f, Metric f, Traversable f, Applicative f, Show (f (f a)), a ~ Double)
       => a                     -- ^ target chi-squared
       -> [Point x a]           -- ^ observations
       -> f (Model x)           -- ^ models
       -> f a                   -- ^ initial weights
       -> [f a]                 -- ^ iterates
maxEnt cAim pts models f = go f
  where
    offTol = 1 -- TODO
    go :: f a -> [f a]
    go f
      | gDiagOff > offTol = error "g not diagonal"
      | otherwise = f' : go f'
      where
        -- determine search subspace
        basis = subspace pts f models  -- <e_i | f_j>
        -- diagonalize m
        m     = basis !*! hessianChiSquared pts models f !*! adjoint basis
        (p, gamma) = eigen3 m

        -- verify simultaneous diagonalization of g
        g     = basis !*! kronecker f !*! adjoint basis
        gDiag = p !*! g !*! adjoint p
        gDiagOff = sum (fmap (sum . fmap (^(2::Int))) gDiag) - norm (diagonal gDiag)

        -- eigenvalues of g
        wG = Debug.Trace.trace (showMatrix basis) $ diagonal gDiag

        -- orthogonalize g
        -- this is where we pass to the eigenbasis
        f' = case length $ filter (nearZero . abs) $ toList wG of
               -- full rank
               0 -> let s = fmap (recip . norm) p -- scale to set g to identity
                        basis' = traceShow ("g", kronecker s !*! gDiag) $ adjoint (kronecker s) !*! basis
                    in goSubspace f (adjoint basis' !*! p) gamma

               -- one dependent pair
               {-
               1 -> let basis' = fmap (^._xy) $ basis ^. _xy
                    in goOrtho f (adjoint basis' !*! p)
                               (fmap (^._xy) $ (^._xy) $ basis !*! m !*! basis)
                               ((^._xy) $ basis !*! gamma)
                               (fmap (^._xy) $ p^._xy)
               -}
               _ -> error "Too many dependent eigenvectors"

    -- | Here we optimize the objective within our orthonormal search subspace
    goSubspace :: forall g. (Foldable g, Metric g, Applicative g, Distributive g)
               => f a        -- ^ current iterate
               -> f (g a)    -- ^ subspace basis
               -> g a        -- ^ gamma, eigenvalues of M
               -> f a
    goSubspace f basis gamma = f'
      where
        -- next iterate
        f'    = f ^+^ basis !* step alpha

        -- limit of l = |df|
        l0    = 0.1 * sum f

        -- various useful quantities
        s0    = entropy f
        s     = adjoint basis !* gradEntropy f
        c0    = chiSquared pts models f
        c     = adjoint basis !* gradChiSquared pts models f

        -- quadratic model of C
        cQ :: g a -> a
        cQ x  = c0 + c `dot` x + sum ((\x g -> g*x*x) <$> x <*> gamma) / 2
        cAim' = max cAim $ c0 - sum ((\x g -> x*x/g) <$> c <*> gamma) / 3

        -- step with Lagrange multiplier alpha
        step :: a -> g a
        step alpha = (\s c g->(alpha * s - c) / (g + alpha))
                     <$> s <*> c <*> gamma

        -- search for alpha to satisfy C == Caim
        alpha = traceShow ("asdf", toDouble c0, cAim') $ findAlpha 1
        findAlpha :: RealFloat a => a -> a
        findAlpha alpha
          | l > l0        = alpha
          | converged     = alpha
          | cV < cAim'    = findAlpha (1.01 * alpha)
          | otherwise     = findAlpha (0.99* alpha)
          where
            x = step alpha
            l = norm x
            cV = cQ x
            converged = (<1e-5) $ abs (cV - cAim') / cAim'
{-# INLINE maxEnt #-}

toDouble = unsafeCoerce :: a -> Double

-- | Component-wise multiplication
cmult :: (Num a, Applicative f) => f a -> f a -> f a
cmult a b = (*) <$> a <*> b
{-# INLINE cmult #-}

cmultOp :: (Num a, Applicative f) => f a -> f (f a) -> f (f a)
cmultOp v a = (*^) <$> v <*> a
{-# INLINE cmultOp #-}

-- | Compute a three-dimensional search subspace
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
{-# INLINE subspace #-}

eigen3 :: M33 Double -> (M33 Double, V3 Double)
eigen3 a =
  let (lambd, p) = LA.eigS $ m33ToHMat a
  in (hmatToM33 p, hvecToV3 lambd)
{-# INLINE eigen3 #-}

m33ToHMat :: HM.Element a => M33 a -> HM.Matrix a
m33ToHMat = HM.fromLists . toList . fmap toList
{-# INLINE m33ToHMat #-}

hmatToM33 :: HM.Element a => HM.Matrix a -> M33 a
hmatToM33 = listToV3 . map listToV3 . HM.toLists
{-# INLINE hmatToM33 #-}

listToV3 :: [a] -> V3 a
listToV3 [x,y,z] = V3 x y z
{-# INLINE listToV3 #-}

hvecToV3 :: Storable a => HV.Vector a -> V3 a
hvecToV3 = listToV3 . HV.toList
{-# INLINE hvecToV3 #-}
