-- | Bunse-Gerstner simultaneous diagonalization

module BunseGurstner where
-- (simultDiag, offAxisTerminate, Diagonalization(..))

import Prelude hiding (sum)
import Linear
import Control.Lens
import Data.Ord
import Data.Complex hiding (conjugate)
import Data.Foldable
import Debug.Trace

data Diagonalization f a = Diag { diagA, diagB, diagQ :: !(f (f a)) }

-- | Simultaneous diagonalization of two normal matricies.
--
-- Requires that,
--     A B = B A
--     A A^H = A^H A
--     B B^H = B^H B
--
-- Method due to,
--
--       A. Bunse-Gerstner, R. Byers, V. Mehrmann. "Numerical Methods for
--       Simultaneous Diagonalization." 1993
simultDiag :: (Conjugate a, RealFloat a, Show a)
           => M33 a
           -> M33 a
           -> [Diagonalization V3 (Complex a)]
simultDiag a b = iterate go (Diag (toComplex a) (toComplex b) eye3)
  where
    toComplex :: (RealFloat a, Functor f) => f (f a) -> f (f (Complex a))
    toComplex = fmap (fmap realToFrac)
    
    go d = foldl' rotPlane d [(ex,ey), (ex,ez), (ey,ez)] 

    rotPlane :: (Conjugate a, RealFloat a, Show a)
             => Diagonalization V3 (Complex a)    -- ^ initial diagonalization
             -> (E V3, E V3)                      -- ^ plane
             -> Diagonalization V3 (Complex a)
    rotPlane (Diag a b q) (i,j) =
        Diag (adjoint r !*! a !*! r)
             (adjoint r !*! b !*! r)
             (q !*! r)
      where
        --(c, s) = rotationAngle a b (i,j)
        (c, s) = bruteForceAngle 100 100 a b (i,j)
        r = rotation i j c s

-- | Terminate the diagonalization when the off-axis norm has pass
-- below epsilon
offAxisTerminate :: (Conjugate a, RealFloat a)
                 => a                                -- ^ epsilon
                 -> [Diagonalization V3 (Complex a)] -- ^ iterates
                 -> [Diagonalization V3 (Complex a)]
offAxisTerminate eps = takeWhile f
  where f d = let off = offAxisNorm (diagA d) + offAxisNorm (diagB d)
                  mag = frobeniusNorm (diagA d) + frobeniusNorm (diagB d)
              in off > eps * mag
 
-- | Weight of components off-axis
offAxisNorm :: (Additive f, Traversable f, Trace f, RealFloat a)
            => f (f (Complex a)) -> a
offAxisNorm = frobeniusNorm . offAxis
  where
    offAxis a = a !-! kronecker (diagonal a)
 
-- | a rotation in ij plane
rotation :: (Conjugate a, RealFloat a, Show a)
         => E V3 -> E V3 -> a -> Complex a -> M33 (Complex a)
rotation _ _ c s | (>1e-5) $ c^2 + (magnitude s)^2 - 1 =
    error "Invalid rotation"
rotation (E ei) (E ej) c s =
    eye3
    + realToFrac (c - 1) *!! (unit ei `outer` unit ei)
    - conjugate s        *!! (unit ei `outer` unit ej)
    +           s        *!! (unit ej `outer` unit ei)
    + realToFrac (c - 1) *!! (unit ej `outer` unit ej)

-- | Pick a rotation
rotationAngle :: (Conjugate a, Ord a, RealFloat a)
              => M33 (Complex a)
              -> M33 (Complex a)
              -> (E V3, E V3)
              -> (a, Complex a)
rotationAngle a b plane =
      minimumBy (comparing $ cond a b plane)
    $ map (\f->f a plane)
    -- $ [goldstineAngle, eberleinAngle]
    $ [eberleinAngle]

bruteForceAngle :: (Conjugate a, Ord a, RealFloat a)
                => Int -> Int
                -> M33 (Complex a)
                -> M33 (Complex a)
                -> (E V3, E V3)
                -> (a, Complex a)
bruteForceAngle nTheta nPhi a b plane =
    minimumBy (comparing $ cond a b plane) $ do
        theta <- range (-pi/4) (pi/4) nTheta
        phi <- range (-pi) pi nPhi
        return (cos theta, cis phi * realToFrac (sin theta))
  where
    range :: RealFrac a => a -> a -> Int -> [a]
    range a b n = [a + d*realToFrac i | i <- [0..n]]
      where d = (b - a) / realToFrac n
    
    
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

-- | Minimizing rotation angle due to Eberlein, 1962
eberleinAngle :: (Conjugate a, Ord a, RealFloat a)
              => M33 (Complex a)
              -> (E V3, E V3)
              -> (a, Complex a)
eberleinAngle a (E i, E j) = (c, realToFrac s)
  where
    -- FIXME? real?
    x = atan2 (realPart $ a^.i.j + a^.j.i) (realPart $ a^.i.i - a^.j.j) / 2
    c = cos x + 1
    s = sin x
  
-- | Minimizing rotation angle due to Goldstine & Horwitz, 1958 (pg. 179)
goldstineAngle :: (Conjugate a, Ord a, RealFloat a)
               => M33 (Complex a)
               -> (E V3, E V3)
               -> (a, Complex a)
goldstineAngle a (E i, E j) = undefined
  where
    -- Hermitian
    a = undefined

frobeniusNorm :: (Foldable f, Functor f, RealFloat a) => f (f (Complex a)) -> a
frobeniusNorm a = sum $ fmap (sum . fmap (\x->(magnitude x)^2)) a
