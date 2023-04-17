-- ----------
-- author:  egor lappo
-- date:    2023-03-11
-- ----------

-- some notes:
-- turtle's `Shell a` is a list monad, so it's necessary to split the Shell actions into little bits
-- if you write the whole script as a single Shell action, it won't work (will go through all possible combinations of all variables bound with <-)
-- wrapping bcftools and beagle into functions is hard, maybe there is a nicer way to do it
-- most of the code is creating/splitting file names, otherwise it is actually really simple
{-# LANGUAGE NamedFieldPuns    #-}
{-# LANGUAGE OverloadedStrings #-}

module Main (main) where

import           Prelude                             hiding (FilePath)
import           Turtle                              hiding (option, options,
                                                      switch)

import           Control.Concurrent.ParallelIO.Local (parallel_, withPool)
import qualified Data.Foldable                       as F
import           Data.List                           (intercalate)
import           Data.Tuple                          (swap)

import           Data.List.Split
import qualified Data.Sequence                       as S
import qualified Data.Text                           as T
import qualified Data.Text.IO                        as T
import           Options.Applicative                 hiding (command)
import           System.Random


data Options = Options {
    skipSplits'     :: !Bool
  , nThreads'       :: !Int
  , nReps'          :: !Int
  , seed'           :: !Int
  , refSize'        :: !Int
  , codisLoci'      :: !FilePath
  , codisChrs'      :: !FilePath
  , sampleIds'      :: !FilePath
  , vcfFolder'      :: !FilePath
  , repsFolder'     :: !FilePath
  , outputFolder'   :: !FilePath
  , beagleJar'      :: !FilePath
  , bref3Jar'       :: !FilePath
  , plinkMapFolder' :: !FilePath
  }

opts :: ParserInfo Options
opts = info (helper <*> options') (fullDesc <> progDesc "run data generation and imputation")
  where
    options' = Options
      <$> switch (long "skip-splits" <> help "Skip writing the replicate splits (only does imputation by beagle in this case)")
      <*> option auto (long "threads" <> short 't' <> metavar "INT" <> help "Number of threads" <> value 6)
      <*> option auto (long "reps" <> short 'n' <> metavar "INT" <> help "Number of replicates" <> value 100)
      <*> option auto (long "seed" <> short 's' <> metavar "INT" <> help "Seed for random number generation" <> value 33)
      <*> option auto (long "refsize" <> short 'r' <> metavar "INT" <> help "Number of samples to be used for training" <> value 1502)
      <*> strOption (long "loci" <> short 'l' <> metavar "FILE" <> help "File containing the list of STR loci" <> value "raw_data/STR_loci.csv")
      <*> strOption (long "chrs" <> short 'c' <> metavar "FILE" <> help "File containing the list of chromosome assignments for STR loci" <> value "raw_data/STR_loci_chr.csv")
      <*> strOption (long "samples" <> short 'i' <> metavar "FILE" <> help "File containing the list of sample ids present in the VCF files" <> value "raw_data/sample_ids.csv")
      <*> strOption (long "vcf" <> short 'v' <> metavar "FILE" <> help "Folder containing the VCF files" <> value "raw_data/vcf/")
      <*> strOption (long "data" <> short 'd' <> metavar "FILE" <> help "Folder to save the replicate splits to" <> value "processed_data/")
      <*> strOption (long "output" <> short 'o' <> metavar "FILE" <> help "Folder to save the results to" <> value "imputation_results/")
      <*> strOption (long "beagle" <> metavar "FILE" <> help "BEAGLE jar file" <> value "beagle/beagle.22Jul22.46e.jar")
      <*> strOption (long "bref" <> metavar "FILE" <> help "bref3 jar file" <> value "beagle/bref3.22Jul22.46e.jar")
      <*> strOption (long "plink" <> metavar "FILE" <> help "Folder containing the plink map files" <> value "beagle/plink/")


-- SCRIPT
main :: IO ()
main = do
  -- read in options, bind to variables
  Options skipSplits threads reps seed' refSize codisLoci codisChrs sampleIds vcfFolder dataFolder outputFolder beagleJar bref3Jar plinkMapFolder <- execParser opts

  -- misc variables
  let
    -- naming helper (extracts locus name from the vcf file name)
    mkName vcf suff = (head . splitOn "_" . basename) vcf <> "_" <> suff
    -- main rng, used to generate seeds for replicates
    rng = mkStdGen seed'
    seeds = randomRs (0, 1000000) rng
    -- folders to contain replicate splits, and
    folders = map ((dataFolder </>) . T.unpack . format ("rep-"%d)) [1..reps]
    -- folders to contain results for replicates
    resultFolders = map ((outputFolder </>) . T.unpack . format ("rep-"%d)) [1..reps]

  -- WRITE THE VCFS FOR REPLICATES
  if skipSplits then putStrLn "skipping making replicate splits" else withPool threads $ \pool -> parallel_ pool $ (flip map) (zip folders seeds) $ \(folder, seed) -> sh $! do
    -- create folders for replicates
    mktree folder
    -- create splits for each replicate
    writeSplitFiles seed refSize sampleIds folder
    -- get the vcf filepaths
    vcf <- find (suffix ".vcf.gz") vcfFolder
    -- write SNPs
    writeSubsetSNPs (folder </> "ref_ids.csv")  vcf (folder </> mkName vcf "ref_snps.vcf.gz")
    writeSubsetSNPs (folder </> "base_ids.csv") vcf (folder </> mkName vcf "base_snps.vcf.gz")
    -- index SNPs
    indexVCF (folder </> mkName vcf "ref_snps.vcf.gz")
    indexVCF (folder </> mkName vcf "base_snps.vcf.gz")
    -- write STRs
    writeSubsetSTRs codisLoci (folder </> "ref_ids.csv")  vcf (folder </> mkName vcf "ref_strs.vcf.gz")
    writeSubsetSTRs codisLoci (folder </> "base_ids.csv") vcf (folder </> mkName vcf "base_strs.vcf.gz")
    -- index STRs
    indexVCF (folder </> mkName vcf "ref_strs.vcf.gz")
    indexVCF (folder </> mkName vcf "base_strs.vcf.gz")
    -- concatenate the SNP and STR files for the reference panel
    concatVCFs (folder </> mkName vcf "ref_snps.vcf.gz") (folder </> mkName vcf "ref_strs.vcf.gz") (folder </> mkName vcf "ref.vcf.gz")
    -- index and sort the concatenated file
    indexVCF (folder </> mkName vcf "ref.vcf.gz")
    sortVCF (folder </> mkName vcf "ref.vcf.gz") (folder </> mkName vcf "ref.vcf.gz")
    -- remove the intermediate reference panel files
    rm (folder </> mkName vcf "ref_snps.vcf.gz")
    rm (folder </> mkName vcf "ref_snps.vcf.gz.csi")
    rm (folder </> mkName vcf "ref_strs.vcf.gz")
    rm (folder </> mkName vcf "ref_strs.vcf.gz.csi")

  -- CLEAN UP THE FILES FOR REPLICATES
  if skipSplits then putStrLn "skipping split cleanup" else withPool threads $ \pool -> parallel_ pool $ (flip map) folders $ \folder' -> sh $! do
    -- merge the database true strs, delete the intermediate files
    concatAllVCFs (folder' </> "*_base_strs.vcf.gz") (folder' </> "base_strs.vcf.gz")
    sortVCF (folder' </> "base_strs.vcf.gz") (folder' </> "base_strs.vcf.gz")
    -- batch remove (can't use turtle's `rm` with a wildcard nicely)
    shell (format ("rm "%fp) (folder' </> "*_base_strs.vcf.gz")) empty
    shell (format ("rm "%fp) (folder' </> "*_base_strs.vcf.gz.csi")) empty

  -- RUN BEAGLE IMPUTATION
  loci <- T.lines <$> T.readFile codisLoci
  chrs <- T.lines <$> T.readFile codisChrs
  -- for each replicate folder
  withPool threads $ \pool -> parallel_ pool $ (flip map) (zip folders resultFolders) $ \(folder, resultFolder) -> do
    -- for each codis locus
    F.forM_ (zip loci chrs) $ \(locus, chr) ->
      let
        snpVCF = folder </> T.unpack locus <> "_base_snps.vcf.gz"
        refVCF = folder </> T.unpack locus <> "_ref.vcf.gz"
        refBref = folder </> T.unpack locus <> "_ref.bref3"
        imputedVCF = folder </> T.unpack locus <> "_base_strs_snps_imputed.vcf.gz"
      in sh $! do
        -- index for beagle
        makeBrefIndex bref3Jar refVCF
        -- run beagle imputation
        runBeagleImpute beagleJar plinkMapFolder refBref chr snpVCF
        -- write only the row for the codis locus to a separate file
        writeSTRs codisLoci imputedVCF (folder </> T.unpack locus <> "_base_strs_imputed.vcf.gz")
        -- index (for concatenation below)
        indexVCF (folder </> T.unpack locus <> "_base_strs_imputed.vcf.gz")

    -- CLEAN UP THE FILES POST-IMPUTATION
    sh $! do
      -- gather across loci and write to a single file
      concatAllVCFs (folder </> "*_base_strs_imputed.vcf.gz") (folder </> "base_strs_imputed.vcf.gz")
      -- remove temp files
      shell (format ("rm "%fp) (folder </> "*_base_strs_imputed.vcf.gz")) empty
      shell (format ("rm "%fp) (folder </> "*.log")) empty
      shell (format ("rm "%fp) (folder </> "*.bref3")) empty

      -- move the results to the result folder (to avoid cluttering the directories)
      -- hopefully, i will never need to look into processed_data...
      mktree resultFolder
      cp (folder </> "base_strs.vcf.gz") (resultFolder </> "base_strs.vcf.gz")
      cp (folder </> "base_strs_imputed.vcf.gz") (resultFolder </> "base_strs_imputed.vcf.gz")
      -- index needed for viewing
      indexVCF (resultFolder </> "base_strs.vcf.gz")
      indexVCF (resultFolder </> "base_strs_imputed.vcf.gz")

      -- make tsvs to never go back to vcfs again
      tableize resultFolder

-- FUNCTIONS

-- BCFTOOLS FUNCTIONS

-- concatenate vcfs
concatVCFs
  :: FilePath -- vcf1
  -> FilePath -- vcf2
  -> FilePath -- destination file
  -> Shell ExitCode
concatVCFs vcf1 vcf2 destFile =
  let
    command = format ("bcftools concat -a -O z -o "%fp%" "%fp%" "%fp) destFile vcf1 vcf2
  in do
    printLine $ format ("Concatenating VCFs: "%w) command
    shell command empty

-- concatenate all vcfs with a given suffix
concatAllVCFs
  :: FilePath -- suffix to select files for concatenation
  -> FilePath -- destination file
  -> Shell ExitCode
concatAllVCFs files destFile = do
  printLine $ format ("Concatenating VCFs: "%w) files
  shell (format ("bcftools concat -a -O z -o "%fp%" "%fp) destFile files) empty

-- index a vcf file
indexVCF :: FilePath -> Shell ExitCode
indexVCF vcfFile =
  let
    command = format ("bcftools index "%fp) vcfFile
  in do
    printLine $ format ("Indexing VCF: "%w) command
    shell command empty

-- sort a vcf file
sortVCF :: FilePath -> FilePath -> Shell ExitCode
sortVCF vcfFile destFile =
  let
    command = format ("bcftools sort -O z -o "%fp%" "%fp) destFile vcfFile
  in do
    printLine $ format ("Sorting VCF: "%w) command
    shell command empty

queryVCF
  :: Text -- query
  -> FilePath -- vcf file
  -> FilePath -- output file
  -> Shell ExitCode
queryVCF query vcfFile outFile =
  let
    command = format ("bcftools query -f "%s%" "%fp%" > "%fp) query vcfFile outFile
  in do
    printLine $ format ("Querying VCF: "%w) command
    shell command empty

-- subset snp data and copy to respective folder
writeSubsetSNPs
  :: FilePath -- sample ids to keep
  -> FilePath -- vcf file to copy from
  -> FilePath -- destination file
  -> Shell ExitCode
writeSubsetSNPs sampleIds vcfFile destFile =
  let
    command = format ("bcftools view -S "%fp%" -m 2 -M 2 --types snps -O z -o"%fp%" "%fp) sampleIds destFile vcfFile
  in do
    printLine $ format ("Writing SNP variants: "%w) command
    shell command empty

-- subset STR data and copy to respective folder
writeSubsetSTRs
  :: FilePath -- file with CODIS loci
  -> FilePath -- sample ids to keep
  -> FilePath -- vcf file to copy from
  -> FilePath -- destination file
  -> Shell ExitCode
writeSubsetSTRs codisLoci sampleIds vcfFile destFile =
  let
    command = format ("bcftools view -S "%fp%" --include ID=@"%fp%" -O z -o"%fp%" "%fp) sampleIds codisLoci destFile vcfFile
  in do
    printLine $ format ("Writing STR variants: "%w) command
    shell command empty

-- subset STR data and copy to respective folder
writeSTRs
  :: FilePath -- file containning codis loci names
  -> FilePath -- vcf file to copy from
  -> FilePath -- destination file
  -> Shell ExitCode
writeSTRs codisLoci vcfFile destFile =
  let
    command = format ("bcftools view --include ID=@"%fp%" -O z -o"%fp%" "%fp) codisLoci destFile vcfFile
  in do
    printLine $ format ("Writing STR variants: "%w) command
    shell command empty

-- turn vcfs into tsvs with bcftools query
tableize :: FilePath -> Shell ExitCode
tableize folder =
  let
    strVCF            = folder </> "base_strs.vcf.gz"
    strImputedVCF     = folder </> "base_strs_imputed.vcf.gz"
    strTable        = folder </> "base_strs.tsv"
    strImputedTable = folder </> "base_strs_imputed.tsv"
  in do
    queryVCF "\"[%ID\\t%SAMPLE\\t%GT\\n]\"" strVCF strTable
    queryVCF "\"[%ID\\t%SAMPLE\\t%GT\\t%DS\\t%AP1\\t%AP2\\t%GP\\n]\"" strImputedVCF strImputedTable

-- BEAGLE FUNCTIONS

makeBrefIndex :: FilePath-> FilePath -> Shell ExitCode
makeBrefIndex bref3Jar ref =
  let
    brefIndexFile = (dropExtension . dropExtension) ref <.> "bref3"
    command = format ("java -Xmx16g -jar "%fp%" "%fp%" > "%fp) bref3Jar ref brefIndexFile
  in do
    printLine $ format ("Making .bref3 index for "%fp) ref
    shell command empty

runBeagleImpute
  :: FilePath -- path to beagle jar
  -> FilePath -- path to plink map folder
  -> FilePath -- reference
  -> Text     -- chromosome (e.g. "5")
  -> FilePath -- SNP file
  -> Shell ExitCode
runBeagleImpute beagleJar plinkMapFolder ref chr snpVCF =
  let
    locus = (head . splitOn "_" . basename) snpVCF
    plinkMap = plinkMapFolder </> ("plink.chr" <> T.unpack chr <> ".GRCh37.map")
    beagleOut = directory snpVCF </> locus <> "_base_strs_snps_imputed"
    command = format ("java -Xmx16g -jar "%fp%" ref="%fp%" gt="%fp%" out="%fp%" map="%fp%" ap=true gp=true impute=true nthreads=8") beagleJar ref snpVCF beagleOut plinkMap
  in do
    printLine $ format ("Running Beagle on "%fp%" for chromosome "%s) snpVCF chr
    shell command empty

-- MISC FUNCTIONS

allExtensions :: FilePath -> FilePath
allExtensions = intercalate "." . snd . splitExtensions

-- convert [Text] to Shell Line
catText :: [Text] -> Shell Line
catText = cat . map (return . unsafeTextToLine)

printLine :: Text -> Shell ()
printLine = echo . unsafeTextToLine

-- given an rng and a destination folder, write out the sample ids chosed for imputation reference and for the database itself
writeSplitFiles
  :: Int -- seed for rng
  -> Int -- reference panel size
  -> FilePath -- file with sampleIds
  -> FilePath -- destination folder
  -> Shell ()
writeSplitFiles seed refSize sampleIds folder = do
    -- read in sample ids
    ((refIds, baseIds), _) <- liftIO (splitList (mkStdGen seed) refSize . T.lines <$> T.readFile sampleIds)
    -- write out sample ids
    output (folder </> "ref_ids.csv") (catText refIds)
    output (folder </> "base_ids.csv") (catText baseIds)

-- sample without replacement
-- https://stackoverflow.com/a/13782206
sample
  :: RandomGen g
  => g        -- the generator to use
  -> Int      -- the size of the sample
  -> [a]      -- the list to sample from
  -> ([a], g) -- the sample and the final generator
sample g size ys = go 0 (len-1) (S.fromList ys) g
  where
    len = length ys
    go n i xs rng
      | n >= size = ((F.toList . S.drop (len - size)) xs, rng)
      | otherwise =
          let
            (j, rng') = uniformR (0, i) rng
            toI  = xs `S.index` j
            toJ  = xs `S.index` i
            nxt = (S.update i toI . S.update j toJ) xs
          in go (n + 1) (i - 1) nxt rng'

-- split a list into two random lists of a given size
-- https://stackoverflow.com/a/13782206
splitList
  :: RandomGen g
  => g              -- the generator to use
  -> Int            -- the size of the first list
  -> [a]            -- the list to split
  -> (([a],[a]), g) -- the two lists and the final generator
splitList g size ys = go 0 (len-1) (S.fromList ys) g
  where
    len = length ys
    go n i xs rng
      | n >= size = ((swap . splitAt (len - size) . F.toList) xs, rng)
      | otherwise =
          let
            (j, rng') = uniformR (0, i) rng
            toI  = xs `S.index` j
            toJ  = xs `S.index` i
            nxt = (S.update i toI . S.update j toJ) xs
          in go (n + 1) (i - 1) nxt rng'
