module Main (main) where

import System.Environment (getArgs)
import System.CPUTime
import Control.DeepSeq (deepseq)
import Text.Printf (printf)
import MyLib
import System.Random

main :: IO ()
main = do

  setStdGen (mkStdGen 42)

  size <- readMatrixSize "input_matrix.txt"
  printf "Matrix size read from file: %d\n" size

  (mat, vec) <- generateMatrix size
  
  -- Measure time for parallel algorithm
  startPar <- getCPUTime
  let solutionPar = solvePar mat vec
  solutionPar `deepseq` return ()
  endPar <- getCPUTime


  -- Measure time for sequential algorithm
  -- startSer <- getCPUTime
  -- let solutionSer = solveSer mat vec
  -- solutionSer `deepseq` return ()
  -- endSeq <- getCPUTime

  let diffPar = fromIntegral (endPar - startPar) / (10^12) :: Double
  printf "\nExecution time for parallel algorithm: %.6f\n" diffPar

  -- let diffSeq = fromIntegral (endSeq - startSeq) / (10^12) :: Double
  -- printf "Execution time for sequential algorithm: %.6f\n" diffSeq

