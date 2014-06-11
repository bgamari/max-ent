import Prelude hiding (sum)                
import Data.Foldable
import Data.Complex
import Control.Applicative
import Linear

import MaxEnt

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
    print $ maxEnt 1 pts models (V4 0.2 0.3 0.2 0.3)
    
