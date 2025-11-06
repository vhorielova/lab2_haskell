module Main (main) where

import System.Environment (getArgs)
import MyLib

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
