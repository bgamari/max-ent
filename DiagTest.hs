import Prelude hiding (sum)                
import Data.Foldable
import Data.Complex
import Control.Applicative
import Linear
import Data.List (intercalate)

import BunseGurstner

main = do
    mapM (\(i,d)->do
      print ("hello", i, offAxisNorm $ diagA d, offAxisNorm $ diagB d)
      putStrLn $ showMatrix $ diagQ d !*! adjoint (diagQ d)
      --putStrLn $ showMatrix $ adjoint (diagQ d) !*! fmap (fmap realToFrac) m !*! diagQ d
      --putStrLn $ showMatrix $ diagQ d !*! fmap (fmap realToFrac) m !*! adjoint (diagQ d)
      putStrLn $ showMatrix $ diagA d
      putStrLn $ showMatrix $ diagB d
      mapRotation (show i) (diagA d) (diagB d)
      )
      $ zip [0..] $ take 40 $ drop 00000
      $ offAxisTerminate 1e-7
      $ simultDiag n m

toComplex :: RealFloat a => M33 a -> M33 (Complex a)
toComplex = fmap (fmap realToFrac)

mapRotation :: (Conjugate a, RealFloat a, Show a)
            => String -> M33 (Complex a) -> M33 (Complex a) -> IO ()
mapRotation t a b = do
    writeRot ("xy-"++t++".txt") (ex,ey)
    writeRot ("xz-"++t++".txt") (ex,ez)
    writeRot ("yz-"++t++".txt") (ey,ez)
    return ()
  where
    writeRot name plane = 
      writeFile name $ unlines
      $ map (\(a,b,c)->intercalate "\t" $ map show [a,b,c]) $ values a b plane
    nTheta = 100
    nPhi = 100

    range :: RealFrac a => a -> a -> Int -> [a]
    range a b n = [a + d*realToFrac i | i <- [0..n]]
      where d = (b - a) / realToFrac n

    --values :: RealFrac a => M33 a -> M33 a -> (E V3, E V3) -> [(a, a, a)]
    values a b plane = do
      theta <- range (-pi/4) (pi/4) nTheta
      phi <- range (-pi) pi nPhi
      let rot = (cos theta, cis phi * realToFrac (sin theta))
      let v = cond a b plane rot
      return (theta, phi, v)
    
showMatrix :: (Functor f, Foldable f, Show a) => f (f a) -> String
showMatrix a = unlines $ toList $ fmap (foldMap (pad 70 . show)) a
  where pad n s = take n $ s ++ repeat ' '

u = axisAngle (V3 0 1 0) (35*pi/180)

vs = fromQuaternion u !*! eye3

n,m :: M33 Double
m = vs !*! kronecker (V3 80 6 1) !*! adjoint vs
n = vs !*! kronecker (V3 7 2 3) !*! adjoint vs

preconditions a b = 
    [ norm (commutator a b)
    , norm $ a !*! adjoint a !-! adjoint a !*! a
    , norm $ b !*! adjoint b !-! adjoint b !*! b
    ]
  where norm = sum . fmap (sum . fmap (^2))
    
commutator a b = a !*! b !-! b !*! a
