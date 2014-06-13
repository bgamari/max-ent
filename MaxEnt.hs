{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ImpredicativeTypes #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell #-}

module MaxEnt where
-- ( Point(..)
-- , Model(..)
-- , mixtureModel
-- , entropy
-- , chiSquared
-- , maxEnt
-- , test
-- ) where

import Prelude hiding (sum, any, minimum)
import Data.Foldable as F
import Data.Traversable as T
import Data.Distributive
import Data.Monoid
import Control.Applicative
import Control.Lens
import Control.Monad.State.Strict
import Linear

import qualified Data.Packed.Matrix as HM
import qualified Data.Packed.Vector as HV
import qualified Numeric.LinearAlgebra.LAPACK as LA
import Foreign.Storable

import Debug.Trace
import Unsafe.Coerce
import Data.List (intercalate)

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

l1Normalize :: (Functor f, Foldable f, RealFrac a)
            => f a -> f a
l1Normalize x = fmap (/ norm) x
  where norm = sum x
{-# INLINE l1Normalize #-}

entropy :: (Foldable f, Metric f, Epsilon a, RealFloat a)
        => f a -> a
entropy = getSum . F.foldMap (\p->Sum $ p * log p) . l1Normalize
{-# INLINE entropy #-}

gradEntropy :: (Foldable f, Functor f, Metric f, Epsilon a, RealFloat a)
            => f a -> f a
gradEntropy = fmap (\p->1 + log p) . l1Normalize
{-# INLINE gradEntropy #-}

hessianEntropy :: (Traversable f, Metric f, Epsilon a, RealFloat a)
               => f a -> f (f a)
hessianEntropy = kronecker . fmap recip . l1Normalize
{-# INLINE hessianEntropy #-}


chiSquared :: (Foldable f, Applicative f, RealFloat a)
           => [Point x a]     -- ^ points
           -> f (Model x)     -- ^ models
           -> f a             -- ^ weights
           -> a
chiSquared pts models f = sum $ fmap g pts
  where
    g (Point x y s) = (mixtureModel models f x - y)^2 / s^2
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


showMatrix :: (Functor f, Foldable f, Functor g, Foldable g, Show a) => f (g a) -> String
showMatrix = (++"\n") . bracket . intercalate ",\n" . toList . fmap row
  where
    pad n s = take n $ s ++ repeat ' '
    bracket c = "[ "++c++" ]"
    row = bracket . intercalate ", " . toList . fmap (pad 30 . show)

tr = Debug.Trace.trace
trm msg = tr msg . tr
trf msg f x = Debug.Trace.trace (show msg++f x) x

data ChopState g a = Chop { _chopAlpha, _chopP :: a
                          , _chopLastX :: g a }

makeLenses ''ChopState

-- | Maximum entropy fitting of distribution
maxEnt :: forall f x a.
          (Distributive f, Metric f, Traversable f, Applicative f, Show (f a), Show (f (f a)), a ~ Double)
       => a                     -- ^ target chi-squared
       -> [Point x a]           -- ^ observations
       -> f (Model x)           -- ^ models
       -> f a                   -- ^ initial weights
       -> [f a]                 -- ^ iterates
maxEnt cAim pts models f = go f
  where
    offTol = 1e-2 -- fraction of frobenius norm of g that can be off-diagonal
    go :: f a -> [f a]
    go f
      | gDiagOff > offTol = error "g not diagonal"
      | otherwise = trm "g" (showMatrix g) $ trm "basis" (showMatrix basis) $ f' : go f'
      where
        -- determine search subspace
        basis = subspace pts f models  -- <e_i | f_j>
        -- diagonalize m
        m     = basis !*! hessianChiSquared pts models f !*! adjoint basis
        (p, gamma) = trf "eig" show $ eigen3 m

        -- verify simultaneous diagonalization of g
        g     = basis !*! kronecker f !*! adjoint basis
        gDiag = adjoint p !*! g !*! p
        gDiagOff = let all = sum $ fmap (sum . fmap (^(2::Int))) gDiag
                       diag = quadrance (diagonal gDiag)
                   in 1 - diag / all

        -- eigenvalues of g
        wG = diagonal gDiag

        -- orthogonalize g
        -- this is where we pass to the eigenbasis
        df = case length $ filter ((< 1e-3) . abs) $ toList wG of
               -- full rank
               0 -> let s = fmap sqrt wG -- scale to set g to identity
                        --p' = kronecker (fmap recip s) !*! p
                        p' = p
                        x = goSubspace f (adjoint basis !*! p') gamma
                    in adjoint basis !* x

               -- one dependent pair
               1 -> let s = fmap sqrt $ wG ^.  _xy
                        --p' = kronecker (fmap recip s) !*! 
                        p' = fmap (^. _xy) (p ^. _xy)
                        basis' = basis ^.  _xy
                        x = goSubspace f (adjoint basis' !*! p')
                            (gamma ^. _xy)
                    in adjoint basis' !* x

               _ -> error "Too many dependent eigenvectors"

        -- next iterate (FIXME: revisit non-negative)
        f'    = fmap (max 1e-8) $ f ^+^ df

    -- | Here we optimize the objective within our orthonormal search subspace
    goSubspace :: forall g. (Foldable g, Metric g, Applicative g, Distributive g, Show (g a))
               => f a        -- ^ current iterate
               -> f (g a)    -- ^ subspace basis
               -> g a        -- ^ gamma, eigenvalues of M
               -> g a
    goSubspace f basis gamma = x
      where
        -- limit of l = |df|
        --l0    = sqrt $ 0.1 * sum f
        l0    = sqrt $ 0.1 * sum (adjoint basis !* f)

        -- various useful quantities
        s     = adjoint basis !* gradEntropy f
        c0    = chiSquared pts models f
        c     = adjoint basis !* gradChiSquared pts models f

        -- @cP p x@ is a quadratic model of C at point x with distance penalty p
        cP :: a -> g a -> a
        cP p x = c0 + c `dot` x + sum ((\x g -> (p+g)*x*x) <$> x <*> gamma) / 2
        cMin = c0 - sum ((\x g -> x*x/g) <$> c <*> gamma) / 2
        cAim'  = traceShow ("cMin = ", cMin) $ max cAim cMin

        -- step with Lagrange multipliers alpha and p
        step :: a -> a -> g a
        step p alpha = (\s c g->(alpha * s - c) / (g + alpha + p))
                       <$> s <*> c <*> gamma

        -- search for alpha to satisfy C == Caim
        alphaMin = case minimum gamma of
                     x | x < 0  -> -x
                     _          -> 0
        alpha0 = alphaMin + 1
        x = evalState alphaChop $ Chop { _chopLastX = step 0 alpha0
                                       , _chopP     = 0
                                       , _chopAlpha = max alpha0 alphaMin
                                       }

        alphaChop :: RealFloat a => State (ChopState g a) (g a)
        alphaChop = do
          p <- use chopP
          x <- uses chopAlpha $ step p
          let cq = cP 0 x
              cp = cP p x
          case () of
            _ | cp > c0       -> tr "cp>c0" $ chopDown False
              | cq < cAim'    -> tr "c<cAim" $ chopUp False
              | norm x > l0   -> tr "too long" $ chopUp False
              | otherwise     -> do chopLastX .= x
                                    chopDown True
          where
            chopDown success = chopAlpha *= 0.93 >> pChop success
            chopUp   success = chopAlpha *= 1.03 >> pChop success

        pChop :: RealFloat a
              => Bool -- ^ alpha chop successful?
              -> State (ChopState g a) (g a)
        pChop success = do
          p <- use chopP
          x <- uses chopAlpha $ step p
          let cq  = cP 0 x
              cp = cP p x
              aChopDone = cAim' <= cq && cq <= cp && cp <= c0
              pChopDone = p == 0
          a <- use chopAlpha
          traceShow (a, p, norm x, cAim', cq, cp, c0, c `dot` x) $ return ()
          case () of
            _ | not aChopDone -> alphaChop
              | not success   -> increaseP >> alphaChop
              -- | not pChopDone -> decreaseP >> alphaChop
              | otherwise     -> use chopLastX
          where
            increaseP = chopP += 0.1 >> chopAlpha .= alpha0
            --decreaseP = chopP -= 0.1 >> chopAlpha .= alpha0
{-# INLINE maxEnt #-}

test :: (Traversable f, Additive f, Applicative f, Metric f, Epsilon a, RealFloat a)
     => [Point x a] -> f (Model x) -> f a -> a
test pts models f = quadrance (gradS ^/ norm gradS ^-^ gradC ^/ norm gradC) / 2
  where
    gradS = gradEntropy f
    gradC = gradChiSquared pts models f
{-# INLINE test #-}

-- | Component-wise multiplication
cmult :: (Num a, Applicative f) => f a -> f a -> f a
cmult a b = (*) <$> a <*> b
{-# INLINE cmult #-}

cmultOp :: (Num a, Applicative f) => f a -> f (f a) -> f (f a)
cmultOp v a = (*^) <$> v <*> a
{-# INLINE cmultOp #-}

-- | Compute a three-dimensional search subspace
subspace :: (Applicative f, Traversable f, Metric f, RealFloat a, Epsilon a)
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

-- | A vector given in the eigenbasis of an operator of type @f (f a)@
newtype Eigenspace f a = Eig (f a)
--eigen3 :: f (f Double) -> (f (Eigenspace f a), f a)

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
