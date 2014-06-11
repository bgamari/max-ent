import Data.Foldable
import Control.Applicative
import Linear

import Test.QuickCheck.Property
import Test.QuickCheck.Modifiers

import MaxEnt
import BunseGurstner

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
    --print $ eigen3 g
    --print $ eigen3 m
    
    --mapM (putStrLn . showMatrix . (\d->diagQ d !*! diagQ d))
    mapM (\d->do print $ offAxisNorm $ diagA d
                 putStrLn $ showMatrix $ diagQ d !*! adjoint (diagQ d)
                 putStrLn $ showMatrix $ adjoint (diagQ d) !*! diagA d !*! diagQ d
                 putStrLn $ showMatrix $ adjoint (diagQ d) !*! diagB d !*! diagQ d
         )
      $ take 100 $ drop 00000
      $ offAxisTerminate 0.1 $ simultDiag n m

showMatrix :: (Functor f, Foldable f, Show a) => f (f a) -> String
showMatrix a = unlines $ toList $ fmap (foldMap (pad 70 . show)) a
  where pad n s = take n $ s ++ repeat ' '

p :: V3 (V3 Double)
p = V3 (V3 1 0 0)
       (V3 0 1 1)
       (V3 0 1 0)

pinv = case inv33 p of Just a -> a

m,n :: V3 (V3 Double)
(m,n) = case commuting3 1 3 2 of
          m:n:_ -> (her m, her n)
   where her a = a + adjoint a

commutator a b = a !*! b !-! b !*! a

-- | Family of commutative matricies
--
-- From,
--
--     E Pereira, C Rosa. "A method to construct sets of commuting
--     matrices"
commuting3 :: Num a => a -> a -> a -> [M33 a]
commuting3 x1 x2 x3 = iterate (^+^ eye3) a
  where
    a = V3 (V3 (x1+x2-x3)  (4*x3-x2)  (2*x2)   )
           (V3 (x2-2*x3)   (x1-x2)     x3      )
           (V3 (-4*x3)     (4*x3)     (x1+2*x3))
           
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
