# CMake script for Bio++ Sequence Library
# Author: Sylvain Gaillard and Julien Dutheil
# Created: 20/08/2009

# File list
SET(CPP_FILES
  Bpp/Seq/DNAToRNA.cpp
  Bpp/Seq/NucleicAcidsReplication.cpp
  Bpp/Seq/Sequence.cpp
  Bpp/Seq/SequenceWithAnnotation.cpp
  Bpp/Seq/SequenceWithAnnotationTools.cpp
  Bpp/Seq/SequenceWithQuality.cpp
  Bpp/Seq/SequenceExceptions.cpp
  Bpp/Seq/SequencePositionIterators.cpp
  Bpp/Seq/SequenceTools.cpp
  Bpp/Seq/SequenceWithQualityTools.cpp
  Bpp/Seq/SequenceWalker.cpp
  Bpp/Seq/Site.cpp
  Bpp/Seq/SiteExceptions.cpp
  Bpp/Seq/SiteTools.cpp
  Bpp/Seq/StringSequenceTools.cpp
  Bpp/Seq/SymbolList.cpp
  Bpp/Seq/SymbolListTools.cpp
  Bpp/Seq/Transliterator.cpp
  Bpp/Seq/DistanceMatrix.cpp
  Bpp/Seq/CodonSiteTools.cpp
  Bpp/Seq/Alphabet/AbstractAlphabet.cpp
  Bpp/Seq/Alphabet/LetterAlphabet.cpp
  Bpp/Seq/Alphabet/AlphabetExceptions.cpp
  Bpp/Seq/Alphabet/AlphabetTools.cpp
  Bpp/Seq/Alphabet/CodonAlphabet.cpp
  Bpp/Seq/Alphabet/DefaultAlphabet.cpp
  Bpp/Seq/Alphabet/CaseMaskedAlphabet.cpp
  Bpp/Seq/Alphabet/DNA.cpp
  Bpp/Seq/Alphabet/ProteicAlphabet.cpp
  Bpp/Seq/Alphabet/RNA.cpp
  Bpp/Seq/Alphabet/WordAlphabet.cpp
  Bpp/Seq/Alphabet/RNY.cpp
  Bpp/Seq/Alphabet/BinaryAlphabet.cpp
  Bpp/Seq/Alphabet/IntegerAlphabet.cpp
  Bpp/Seq/Alphabet/NumericAlphabet.cpp
  Bpp/Seq/GeneticCode/GeneticCode.cpp
  Bpp/Seq/GeneticCode/StandardGeneticCode.cpp
  Bpp/Seq/GeneticCode/VertebrateMitochondrialGeneticCode.cpp
  Bpp/Seq/GeneticCode/InvertebrateMitochondrialGeneticCode.cpp
  Bpp/Seq/GeneticCode/EchinodermMitochondrialGeneticCode.cpp
  Bpp/Seq/GeneticCode/YeastMitochondrialGeneticCode.cpp
  Bpp/Seq/GeneticCode/AscidianMitochondrialGeneticCode.cpp
  Bpp/Seq/GeneticCode/MoldMitochondrialGeneticCode.cpp
  Bpp/Seq/AlphabetIndex/GranthamAAChemicalDistance.cpp
  Bpp/Seq/AlphabetIndex/MiyataAAChemicalDistance.cpp
  Bpp/Seq/AlphabetIndex/BLOSUM50.cpp
  Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.cpp
  Bpp/Seq/AlphabetIndex/AAIndex1Entry.cpp
  Bpp/Seq/AlphabetIndex/AAIndex2Entry.cpp
  Bpp/Seq/AlphabetIndex/SimpleScore.cpp
  Bpp/Seq/AlphabetIndex/__GranthamMatrixCode
  Bpp/Seq/AlphabetIndex/__MiyataMatrixCode
  Bpp/Seq/AlphabetIndex/__BLOSUM50MatrixCode
  Bpp/Seq/Container/AbstractSequenceContainer.cpp
  Bpp/Seq/Container/AlignedSequenceContainer.cpp
  Bpp/Seq/Container/MapSequenceContainer.cpp
  Bpp/Seq/Container/SequenceContainerTools.cpp
  Bpp/Seq/Container/SiteContainerExceptions.cpp
  Bpp/Seq/Container/SiteContainerTools.cpp
  Bpp/Seq/Container/VectorSequenceContainer.cpp
  Bpp/Seq/Container/VectorSiteContainer.cpp
  Bpp/Seq/Container/CompressedVectorSiteContainer.cpp
  Bpp/Seq/Container/SiteContainerIterator.cpp
  Bpp/Seq/Container/SequenceContainerIterator.cpp
  Bpp/Seq/Io/Clustal.cpp
  Bpp/Seq/Io/Dcse.cpp
  Bpp/Seq/Io/Fasta.cpp
  Bpp/Seq/Io/Mase.cpp
  Bpp/Seq/Io/MaseTools.cpp
  Bpp/Seq/Io/Phylip.cpp
  Bpp/Seq/Io/IoSequenceFactory.cpp
  Bpp/Seq/Io/GenBank.cpp
  Bpp/Seq/Io/PhredPhd.cpp
  Bpp/Seq/Io/PhredPoly.cpp
  Bpp/Seq/Io/NexusIoSequence.cpp
  Bpp/Seq/Io/NexusTools.cpp
  Bpp/Seq/Io/Stockholm.cpp
  Bpp/Seq/Io/StreamSequenceIterator.cpp
  Bpp/Seq/Io/BppOSequenceReaderFormat.cpp
  Bpp/Seq/Io/BppOSequenceWriterFormat.cpp
  Bpp/Seq/Io/BppOAlignmentReaderFormat.cpp
  Bpp/Seq/Io/BppOAlignmentWriterFormat.cpp
  Bpp/Seq/Io/BppOSequenceStreamReaderFormat.cpp
  Bpp/Seq/Io/BppOAlphabetIndex1Format.cpp
  Bpp/Seq/Io/BppOAlphabetIndex2Format.cpp
  Bpp/Seq/App/SequenceApplicationTools.cpp
  )

SET(H_FILES
  Bpp/Seq/CodonSiteTools.h
  Bpp/Seq/DNAToRNA.h
  Bpp/Seq/NucleicAcidsReplication.h
  Bpp/Seq/SequenceExceptions.h
  Bpp/Seq/SequencePositionIterators.h
  Bpp/Seq/Sequence.h
  Bpp/Seq/SequenceIterator.h
  Bpp/Seq/SequenceWithAnnotation.h
  Bpp/Seq/SequenceWithAnnotationTools.h
  Bpp/Seq/SequenceWithQuality.h
  Bpp/Seq/SequenceTools.h
  Bpp/Seq/SequenceWithQualityTools.h
  Bpp/Seq/SequenceWalker.h
  Bpp/Seq/SiteExceptions.h
  Bpp/Seq/Site.h
  Bpp/Seq/SiteIterator.h
  Bpp/Seq/SiteTools.h
  Bpp/Seq/StringSequenceTools.h
  Bpp/Seq/SymbolList.h
  Bpp/Seq/SymbolListTools.h
  Bpp/Seq/Transliterator.h
  Bpp/Seq/DistanceMatrix.h
  Bpp/Seq/Alphabet/AbstractAlphabet.h
  Bpp/Seq/Alphabet/AlphabetExceptions.h
  Bpp/Seq/Alphabet/Alphabet.h
  Bpp/Seq/Alphabet/AlphabetState.h
  Bpp/Seq/Alphabet/AlphabetTools.h
  Bpp/Seq/Alphabet/CodonAlphabet.h
  Bpp/Seq/Alphabet/CaseMaskedAlphabet.h
  Bpp/Seq/Alphabet/DefaultAlphabet.h
  Bpp/Seq/Alphabet/DNA.h
  Bpp/Seq/Alphabet/LetterAlphabet.h
  Bpp/Seq/Alphabet/NucleicAlphabet.h
  Bpp/Seq/Alphabet/NucleicAlphabetState.h
  Bpp/Seq/Alphabet/ProteicAlphabet.h
  Bpp/Seq/Alphabet/ProteicAlphabetState.h
  Bpp/Seq/Alphabet/RNA.h
  Bpp/Seq/Alphabet/WordAlphabet.h
  Bpp/Seq/Alphabet/RNY.h
  Bpp/Seq/Alphabet/BinaryAlphabet.h
  Bpp/Seq/Alphabet/IntegerAlphabet.h
  Bpp/Seq/Alphabet/AlphabetNumericState.h
  Bpp/Seq/Alphabet/NumericAlphabet.h
  Bpp/Seq/GeneticCode/GeneticCode.h
  Bpp/Seq/GeneticCode/StandardGeneticCode.h
  Bpp/Seq/GeneticCode/VertebrateMitochondrialGeneticCode.h
  Bpp/Seq/GeneticCode/InvertebrateMitochondrialGeneticCode.h
  Bpp/Seq/GeneticCode/EchinodermMitochondrialGeneticCode.h
  Bpp/Seq/GeneticCode/YeastMitochondrialGeneticCode.h
  Bpp/Seq/Io/AbstractISequence.h
  Bpp/Seq/Io/AbstractIAlignment.h
  Bpp/Seq/Io/AbstractOSequence.h
  Bpp/Seq/Io/AbstractOAlignment.h
  Bpp/Seq/Io/Clustal.h
  Bpp/Seq/Io/Dcse.h
  Bpp/Seq/Io/Fasta.h
  Bpp/Seq/Io/IOSequence.h
  Bpp/Seq/Io/ISequence.h
  Bpp/Seq/Io/ISequenceStream.h
  Bpp/Seq/Io/IoSequenceStream.h
  Bpp/Seq/Io/Mase.h
  Bpp/Seq/Io/MaseTools.h
  Bpp/Seq/Io/OSequence.h
  Bpp/Seq/Io/OSequenceStream.h
  Bpp/Seq/Io/Phylip.h
  Bpp/Seq/Io/GenBank.h
  Bpp/Seq/Io/PhredPhd.h
  Bpp/Seq/Io/PhredPoly.h
  Bpp/Seq/Io/NexusIOSequence.h
  Bpp/Seq/Io/NexusTools.h
  Bpp/Seq/Io/Stockholm.h
  Bpp/Seq/Io/IoSequenceFactory.h
  Bpp/Seq/Io/SequenceFileIndex.h
  Bpp/Seq/Io/StreamSequenceIterator.h
  Bpp/Seq/Io/BppOSequenceReaderFormat.h
  Bpp/Seq/Io/BppOSequenceWriterFormat.h
  Bpp/Seq/Io/BppOAlignmentReaderFormat.h
  Bpp/Seq/Io/BppOAlignmentWriterFormat.h
  Bpp/Seq/Io/BppOSequenceStreamReaderFormat.h
  Bpp/Seq/Io/BppOAlphabetIndex1Format.h
  Bpp/Seq/Io/BppOAlphabetIndex2Format.h
  Bpp/Seq/Container/VectorSequenceContainer.h
  Bpp/Seq/Container/VectorSiteContainer.h
  Bpp/Seq/Container/CompressedVectorSiteContainer.h
  Bpp/Seq/Container/MapSequenceContainer.h
  Bpp/Seq/Container/OrderedSequenceContainer.h
  Bpp/Seq/Container/AbstractSequenceContainer.h
  Bpp/Seq/Container/AlignedSequenceContainer.h 
  Bpp/Seq/Container/SequenceContainerExceptions.h
  Bpp/Seq/Container/SequenceContainer.h
  Bpp/Seq/Container/SequenceContainerTools.h
  Bpp/Seq/Container/SiteContainerExceptions.h
  Bpp/Seq/Container/SiteContainer.h
  Bpp/Seq/Container/SiteContainerTools.h
  Bpp/Seq/Container/SiteContainerIterator.h
  Bpp/Seq/AlphabetIndex/AlphabetIndex1.h
  Bpp/Seq/AlphabetIndex/AlphabetIndex2.h
  Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h
  Bpp/Seq/AlphabetIndex/GranthamAAChemicalDistance.h
  Bpp/Seq/AlphabetIndex/GranthamAAPolarityIndex.h
  Bpp/Seq/AlphabetIndex/GranthamAAVolumeIndex.h
  Bpp/Seq/AlphabetIndex/MiyataAAChemicalDistance.h
  Bpp/Seq/AlphabetIndex/BLOSUM50.h
  Bpp/Seq/AlphabetIndex/SimpleIndexDistance.h
  Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.h
  Bpp/Seq/AlphabetIndex/AAIndex1Entry.h
  Bpp/Seq/AlphabetIndex/AAIndex2Entry.h
  Bpp/Seq/AlphabetIndex/KleinAANetChargeIndex.h
  Bpp/Seq/AlphabetIndex/AAChargeIndex.h
  Bpp/Seq/AlphabetIndex/SimpleScore.h
  Bpp/Seq/AlphabetIndex/AAChouFasmanAHelixIndex.h
  Bpp/Seq/AlphabetIndex/AAChouFasmanBSheetIndex.h
  Bpp/Seq/AlphabetIndex/AAChouFasmanTurnIndex.h
  Bpp/Seq/AlphabetIndex/AAVolumeIndex.h
  Bpp/Seq/AlphabetIndex/AASurfaceIndex.h
  Bpp/Seq/AlphabetIndex/AAMassIndex.h
  Bpp/Seq/AlphabetIndex/AAChenGuHuangHydrophobicityIndex.h
  Bpp/Seq/AlphabetIndex/AASEASup30Index.h
  Bpp/Seq/AlphabetIndex/AASEAInf10Index.h
  Bpp/Seq/AlphabetIndex/AASEA1030Index.h
  Bpp/Seq/App/SequenceApplicationTools.h
  )

# Build the static lib
ADD_LIBRARY(bppseq-static STATIC ${CPP_FILES})
SET_TARGET_PROPERTIES(bppseq-static
  PROPERTIES OUTPUT_NAME bpp-seq
  CLEAN_DIRECT_OUTPUT 1
  )
TARGET_LINK_LIBRARIES(bppseq-static ${LIBS})

# Build the shared lib
ADD_LIBRARY(bppseq-shared SHARED ${CPP_FILES})
SET_TARGET_PROPERTIES(bppseq-shared
  PROPERTIES OUTPUT_NAME bpp-seq
  CLEAN_DIRECT_OUTPUT 1
  MACOSX_RPATH 1
  VERSION ${BPPSEQ_VERSION}
  SOVERSION ${BPPSEQ_VERSION_MAJOR}
  )
TARGET_LINK_LIBRARIES(bppseq-shared ${LIBS})

# Install libs
INSTALL(TARGETS bppseq-static bppseq-shared DESTINATION lib${LIB_SUFFIX})

# Install headers
INSTALL(DIRECTORY Bpp/ DESTINATION include/Bpp FILES_MATCHING PATTERN "*.h")

# Generate generic include files (.all)
INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_SOURCE_DIR}/genIncludes.sh ${CMAKE_PREFIX_PATH}/include/Bpp)")

