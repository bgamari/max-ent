import Prelude hiding (sum)                
import Data.Foldable
import Data.Complex
import Control.Applicative
import Linear
import Data.List (intercalate)

import BunseGurstner

main = do
    mapM (\d->do print "hello"
                 print $ offAxisNorm $ diagA d
                 putStrLn $ showMatrix $ diagQ d !*! adjoint (diagQ d)
                 --putStrLn $ showMatrix $ adjoint (diagQ d) !*! fmap (fmap realToFrac) m !*! diagQ d
                 --putStrLn $ showMatrix $ diagQ d !*! fmap (fmap realToFrac) m !*! adjoint (diagQ d)
                 putStrLn $ showMatrix $ diagA d
                 putStrLn $ showMatrix $ diagB d
         )
      $ take 20 $ drop 00000
      $ offAxisTerminate 0.01 $ simultDiag n m

toComplex = fmap (fmap realToFrac)

mapRotation = do
    writeRot "xy.txt" (ex,ey)
    writeRot "xz.txt" (ex,ez)
    writeRot "yz.txt" (ey,ez)
  where
    writeRot name plane = 
      writeFile name $ unlines
      $ map (\(a,b,c)->intercalate "\t" $ map show [a,b,c]) $ values plane
    nTheta = 100
    nPhi = 100

    range :: RealFrac a => a -> a -> Int -> [a]
    range a b n = [a + d*realToFrac i | i <- [0..n]]
      where d = (b - a) / realToFrac n

    values :: (E V3, E V3) -> [(Double, Double, Double)]
    values plane = do
      theta <- range (-pi/4) (pi/4) nTheta
      phi <- range (-pi) pi nPhi
      let rot = (cos theta, cis phi * realToFrac (sin theta))
      let v = cond (toComplex n) (toComplex m) plane rot
      return (theta, phi, v)
    
showMatrix :: (Functor f, Foldable f, Show a) => f (f a) -> String
showMatrix a = unlines $ toList $ fmap (foldMap (pad 70 . show)) a
  where pad n s = take n $ s ++ repeat ' '

p :: V3 (V3 Double)
p = V3 (V3 1 1 0)
       (V3 1 3 2)
       (V3 0 2 1)

pinv = case inv33 p of Just a -> a

m,n :: V3 (V3 Double)
(m,n) = case commuting3 2 1 1 of
          m:n:_ -> (herm m, herm n)
   where
     herm a = (a + adjoint a) !!* scale
       where scale = recip $ sum $ fmap (sum . (^2)) a

preconditions a b = 
    [ norm (commutator a b)
    , norm $ a !*! adjoint a !-! adjoint a !*! a
    , norm $ b !*! adjoint b !-! adjoint b !*! b
    ]
  where norm = sum . fmap (sum . fmap (^2))
    
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
           
