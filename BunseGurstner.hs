{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- | Bunse-Gerstner simultaneous diagonalization
--
-- Method due to,
--
--       A. Bunse-Gerstner, R. Byers, V. Mehrmann. "Numerical Methods for
--       Simultaneous Diagonalization." 1993

module BunseGurstner
  ( simultDiag
  , offAxisTerminate, offAxisNorm
  , Diagonalization(..)
  , frobeniusNorm
  , cond
  ) where

import Prelude hiding (sum, concatMap)
import Linear
import Control.Applicative
import Control.Lens
import Data.Ord
import Data.Complex hiding (conjugate)
import Data.Foldable

import System.IO.Unsafe
import Debug.Trace
import Unsafe.Coerce

data Diagonalization f a = Diag { diagA, diagB, diagQ :: !(f (f a)) }

-- | Simultaneous diagonalization of two normal matricies.
--
-- Requires that,
--     A B   = B A
--     A A^H = A^H A
--     B B^H = B^H B
--
simultDiag :: (Conjugate a, RealFloat a, Epsilon a, Show a)
           => M33 a
           -> M33 a
           -> [Diagonalization V3 (Complex a)]
simultDiag a b = iterate go (Diag (toComplex a) (toComplex b) eye3)
  where
    toComplex :: (RealFloat a, Functor f) => f (f a) -> f (f (Complex a))
    toComplex = fmap (fmap realToFrac)

    go d = foldl' rotPlane d $ [(ex,ey), (ex,ez), (ey,ez)]

    rotPlane :: (Conjugate a, RealFloat a, Epsilon a, Show a)
             => Diagonalization V3 (Complex a)    -- ^ initial diagonalization
             -> (E V3, E V3)                      -- ^ plane
             -> Diagonalization V3 (Complex a)
    rotPlane (Diag a b q) (i,j) =
        Diag (adjoint r !*! a !*! r)
             (adjoint r !*! b !*! r)
             (q !*! r)
      where
        (c, s) = rotationAngle a b (i,j)
        --(c, s) = bruteForceAngle 1000 1000 a b (i,j)
        --(c, s) = deepenAngle 100 100 a b (i,j)
        --(c, s) = optimizeAngle a b (i,j)
        r = rotation i j c s
{-# INLINE simultDiag #-}

-- | Terminate the diagonalization when the off-axis norm has pass
-- below epsilon
offAxisTerminate :: (Conjugate a, RealFloat a)
                 => a                                -- ^ epsilon
                 -> [Diagonalization V3 (Complex a)] -- ^ iterates
                 -> [Diagonalization V3 (Complex a)]
offAxisTerminate eps (d:rest)
  | off < eps * mag = d : take 1 rest
  | otherwise       = d : offAxisTerminate eps rest
  where
    off = offAxisNorm (diagA d) + offAxisNorm (diagB d)
    mag = frobeniusNorm (diagA d) + frobeniusNorm (diagB d)
{-# INLINE offAxisTerminate #-}

-- | Weight of components off-axis
offAxisNorm :: (Additive f, Traversable f, Trace f, RealFloat a)
            => f (f (Complex a)) -> a
offAxisNorm = frobeniusNorm . offAxis
  where
    offAxis a = a !-! kronecker (diagonal a)
{-# INLINE offAxisNorm #-}

-- | a rotation in ij plane
rotation :: (Conjugate a, RealFloat a, Traversable f, Applicative f)
         => E f -> E f -> a -> Complex a -> f (f (Complex a))
rotation _ _ c s | (>1e-5) $ c^2 + (magnitude s)^2 - 1 =
    error "Invalid rotation"
rotation (E ei) (E ej) c s =
      ei . ei %~ (+ realToFrac (c - 1))
    $ ei . ej .~ (- conjugate s)
    $ ej . ei .~ s
    $ ej . ej %~ (+ realToFrac (c - 1))
    $ kronecker (pure 1)
{-# INLINE rotation #-}

bruteForceAngle :: (Conjugate a, Ord a, RealFloat a)
                => Int -> Int
                -> M33 (Complex a)
                -> M33 (Complex a)
                -> (E V3, E V3)
                -> (a, Complex a)
bruteForceAngle nTheta nPhi a b plane =
    minimumBy (comparing $ cond a b plane) $ do
        theta <- linSpace (-pi/4) (pi/4) nTheta
        phi <- linSpace (-pi) pi nPhi
        return (cos theta, cis phi * realToFrac (sin theta))
{-# INLINE bruteForceAngle #-}

linSpace :: RealFrac a => a -> a -> Int -> [a]
linSpace a b n = [a + d*realToFrac i | i <- [0..n]]
  where d = (b - a) / realToFrac n
{-# INLINE linSpace #-}
   
trf :: String -> (a -> String) -> a -> a
trf msg f x = Debug.Trace.trace (msg++f x) x
asT = unsafeCoerce :: a -> (Double, Double)

deepenAngle :: forall a. (Conjugate a, Ord a, RealFloat a)
            => Int -> Int
            -> M33 (Complex a)
            -> M33 (Complex a)
            -> (E V3, E V3)
            -> (a, Complex a)
deepenAngle nTheta nPhi a b plane = toSinCos $ trf "sc = " (show . asT)
                                  $ head $ drop 8 $ go 1 (0,0)
  where
    toSinCos :: (a, a) -> (a, Complex a)
    toSinCos (th,ph) = (cos th, cis ph * realToFrac (sin th))

    go :: Int -> (a, a) -> [(a, a)]
    go n (cTheta,cPhi) =
      let m :: (a,a)
          m = minimumBy (comparing $ cond a b plane . toSinCos) $ do
            theta <- linSpace (cTheta - pi/4/n') (cTheta + pi/4/n') nTheta
            phi <- linSpace (cPhi - pi/n') (cPhi + pi/n') nPhi
            return (theta, phi)
          n' = realToFrac n
      in m : go (2*n) m
{-# INLINE deepenAngle #-}

{-
optimizeAngle :: (Conjugate a, Ord a, RealFloat a)
              => M33 (Complex a)
              -> M33 (Complex a)
              -> (E V3, E V3)
              -> (a, Complex a)
optimizeAngle a b plane =
    head $ drop 50 $ gradientDescent f (V2 0 0)
  where
    f (V2 th phi) = cond (fmap (fmap realToFrac) a)
                         (fmap (fmap realToFrac) b) plane
                         (cos th, cis phi * realToFrac (sin th))

gradDescent :: (f a -> a) -> f a -> [f a]
gradDescent f x = go 1 x (finiteDiff 1 x)
  where
    go h x diffLast =

    grad :: V2 a -> State (f a) a
    grad x = do
        h <- get
        kronecker h
-}

-- | Looking for minimizer of this
cond :: (Conjugate a, Ord a, RealFloat a)
     => M33 (Complex a)  -- ^ first matrix
     -> M33 (Complex a)  -- ^ second matrix
     -> (E V3, E V3)     -- ^ plane ij
     -> (a, Complex a)   -- ^ rotation parameters
     -> a                -- ^ f_{ij}(c, s)
cond a b (E i, E j) (c, s) =
    magnitude $ norm
    $ m !* V3 ((realToFrac c)^2)
              (sqrt 2 * realToFrac c * s)
              (s^2)
  where
    aii = a ^. i . i
    aij = a ^. i . j
    ajj = a ^. j . j
    aji = a ^. j . i

    bii = b ^. i . i
    bij = b ^. i . j
    bjj = b ^. j . j
    bji = b ^. j . i

    avg i j = (i - j) / sqrt 2
    conj = conjugate
    m = V4 (V3 (conj aij) (avg (conj aii) (conj ajj)) (-conj aji))
           (V3 (     aji) (avg       aii        ajj ) (-     aij))
           (V3 (conj bij) (avg (conj bii) (conj bjj)) (-conj bji))
           (V3 (     bji) (avg       bii        bjj ) (-     bij))
{-# INLINE cond #-}

-- | Pick a rotation
rotationAngle :: (Conjugate a, Ord a, RealFloat a, Epsilon a)
              => M33 (Complex a)
              -> M33 (Complex a)
              -> (E V3, E V3)
              -> (a, Complex a)
rotationAngle a b plane =
      minimumBy (comparing $ cond a b plane)
    $ concatMap (\f->[f a plane, f b plane])
    -- $ [goldstineAngle, eberleinAngle]
    $ [eberleinAngle]
{-# INLINE rotationAngle #-}

-- | Minimizing rotation angle due to Eberlein, 1962
eberleinAngle :: (Conjugate a, Ord a, RealFloat a)
              => M33 (Complex a)
              -> (E V3, E V3)
              -> (a, Complex a)
eberleinAngle a (E i, E j) = (c, realToFrac s)
  where
    -- FIXME? real?
    x = atan2 (realPart $ a^.i.j + a^.j.i) (realPart $ a^.i.i - a^.j.j) / 2
    c = cos x
    s = sin x
{-# INLINE eberleinAngle #-}

-- | Minimizing rotation angle due to Goldstine & Horwitz, 1958 (pg. 179)
goldstineAngle :: (Conjugate a, Ord a, RealFloat a, Epsilon a)
               => M33 (Complex a)
               -> (E V3, E V3)
               -> (a, Complex a)
goldstineAngle a (E i, E j) = (cos phi, cis alpha * realToFrac (sin phi))
  where
    -- Hermitian
    b = (a !+! adjoint a) !!* recip 2
    r i j = magnitude (b ^. i . j)
    beta i j = phase (b ^. i . j)
    u = r i i - r j j
    -- Skew Hermitian
    c = (a !-! adjoint a) !!* recip 2
    s i j = magnitude (c ^. i . j)
    gamma i j = phase (c ^. i . j)
    v = s i i - s j j
    --
    t = u :+ v
    kappa = u^2 + v^2 - 4*((r i j)^2 * (cos $ beta i j - alpha)^2 + (s i j)^2 * (sin $ gamma i j - alpha)^2)
    lambd = 4 * (r i j * u * cos (beta i j - alpha) + s i j * v * sin (gamma i j - alpha))
    --
    phi = acos (signum kappa) / 4
    alpha = let v = case beta i j - gamma i j of
                      x | nearZero x  -> pi / 2
                        | otherwise   -> x
                num = (r i j + s i j)^2 * sin v
                denom = (r i j + s i j)^2 * (sin v)^2 + (r i j - s i j)^2 * (cos v)^2
            in (asin (-num / sqrt denom) - beta i j - gamma i j) / 2
{-# INLINE goldstineAngle #-}

frobeniusNorm :: (Foldable f, Functor f, RealFloat a) => f (f (Complex a)) -> a
frobeniusNorm a = sum $ fmap (sum . fmap (\x->(magnitude x)^2)) a
{-# INLINE frobeniusNorm #-}
