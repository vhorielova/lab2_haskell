module MyLib where

import Control.Parallel.Strategies    
import Data.List (foldl')
import System.Random (randomRIO)

type Matrix = [[Double]]
type Vector = [Double]

gaussianEliminationPar :: Matrix -> Vector -> (Matrix, Vector)
gaussianEliminationPar mat vec = go 0 mat vec
  where
    n = length mat
    go k a b
      | k >= n = (a, b)
      | otherwise =
          let 
              pivotRow = k + maxIndex (drop k (map (abs . (!! k)) a))
              (a', b') = if pivotRow /= k
                         then (swapRows k pivotRow a, swapElems k pivotRow b)
                         else (a, b)
              
              pivot = (a' !! k) !! k
              _ = if abs pivot < 1e-10
                  then error "Matrix is singular"
                  else ()

              indices = [k+1..n-1]
              results = parMap rdeepseq 
                        (\i -> eliminateRow i k pivot (a' !! i) (a' !! k) (b' !! i) (b' !! k))
                        indices
              
              (newRows, newB) = unzip results
              a'' = take k a' ++ [a' !! k] ++ newRows ++ drop n a'
              b'' = take k b' ++ [b' !! k] ++ newB ++ drop n b'
              
          in go (k + 1) a'' b''

gaussianEliminationSer :: Matrix -> Vector -> (Matrix, Vector)
gaussianEliminationSer mat vec = go 0 mat vec
  where
    n = length mat
    go k a b
      | k >= n = (a, b)
      | otherwise =
          let 
              pivotRow = k + maxIndex (drop k (map (abs . (!! k)) a))
              (a', b') = if pivotRow /= k
                         then (swapRows k pivotRow a, swapElems k pivotRow b)
                         else (a, b)
              
              pivot = (a' !! k) !! k
              _ = if abs pivot < 1e-10
                  then error "Matrix is singular"
                  else ()

              indices = [k+1..n-1]
              results = 
                  [ eliminateRow i k pivot (a' !! i) (a' !! k) (b' !! i) (b' !! k)
                  | i <- indices ]
              
              (newRows, newB) = unzip results
              a'' = take (k + 1) a' ++ newRows
              b'' = take (k + 1) b' ++ newB
          in go (k + 1) a'' b''


eliminateRow :: Int -> Int -> Double -> Vector -> Vector -> Double -> Double -> (Vector, Double)
eliminateRow i k pivot rowI pivotRow bi bk =
  let factor = (rowI !! k) / pivot
      newRow = zipWith (\x y -> x - factor * y) rowI pivotRow
      newB = bi - factor * bk
  in (newRow, newB)

backSubstitution :: Matrix -> Vector -> Vector
backSubstitution mat vec = go (n - 1) (replicate n 0)
  where
    n = length mat
    
    go (-1) sol = sol
    go i sol =
      let row = mat !! i
          s = sum $ zipWith (*) (drop (i + 1) row) (drop (i + 1) sol)
          xi = (vec !! i - s) / (row !! i)
          sol' = take i sol ++ [xi] ++ drop (i + 1) sol
      in go (i - 1) sol'

solvePar :: Matrix -> Vector -> Vector
solvePar mat vec = 
  let (mat', vec') = gaussianEliminationPar mat vec
  in backSubstitution mat' vec'

solveSer :: Matrix -> Vector -> Vector
solveSer mat vec =
    let (mat', vec') = gaussianEliminationSer mat vec
    in backSubstitution mat' vec'

maxIndex :: [Double] -> Int
maxIndex xs = snd $ foldl' (\(m, mi) (x, i) -> if x > m then (x, i) else (m, mi)) 
                           (head xs, 0) 
                           (zip xs [0..])

swapRows :: Int -> Int -> Matrix -> Matrix
swapRows i j mat =
  let ri = mat !! i
      rj = mat !! j
  in take i mat ++ [rj] ++ take (j - i - 1) (drop (i + 1) mat) ++ 
     [ri] ++ drop (j + 1) mat

swapElems :: Int -> Int -> [a] -> [a]
swapElems i j xs =
  let ei = xs !! i
      ej = xs !! j
  in take i xs ++ [ej] ++ take (j - i - 1) (drop (i + 1) xs) ++ 
     [ei] ++ drop (j + 1) xs

readMatrixSize :: FilePath -> IO Int
readMatrixSize path = do
  content <- readFile path
  let n = read (head (lines content)) :: Int
  return n

generateMatrix :: Int -> IO (Matrix, Vector)
generateMatrix n = do
    a <- mapM (\_ -> generateRow n) [1..n]
    b <- mapM (\_ -> randomRIO (-100, 100)) [1..n] 
    return (a, b)


generateRow :: Int -> IO [Double]
generateRow n = mapM (\_ -> randomRIO (-100, 100)) [1..n]

