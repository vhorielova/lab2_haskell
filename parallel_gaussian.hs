import Control.Parallel.Strategies
import Control.DeepSeq
import Data.List (foldl')
import System.Environment (getArgs)

type Matrix = [[Double]]
type Vector = [Double]

gaussianElimination :: Matrix -> Vector -> (Matrix, Vector)
gaussianElimination mat vec = go 0 mat vec
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

solve :: Matrix -> Vector -> Vector
solve mat vec = 
  let (mat', vec') = gaussianElimination mat vec
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

readMatrix :: FilePath -> IO (Matrix, Vector)
readMatrix path = do
  content <- readFile path
  let (nStr:rows) = lines content
      n = read nStr :: Int
      parsed = map (map read . words) rows :: [[Double]]
      a = map init parsed
      b = map last parsed
  return (a, b)

main :: IO ()
main = do
  args <- getArgs
  let size = if null args then 4 else read (head args) :: Int

  (mat, vec) <- readMatrix "input_matrix.txt"
  
  putStrLn "Solving system Ax = b"
  putStrLn "Matrix A:"
  mapM_ print mat
  putStrLn "\nVector b:"
  print vec
  
  let solution = solve mat vec
  putStrLn "\nSolution x:"
  print solution
  
  let result = map (\row -> sum $ zipWith (*) row solution) mat
  putStrLn "\nVerification (Ax):"
  print result
  putStrLn "\nExpected (b):"
  print vec
