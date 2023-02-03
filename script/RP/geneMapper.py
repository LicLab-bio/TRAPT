#130428

#ROSE_geneMapper.py

import os

from tqdm import tqdm

#main method wrapped script to take the enhancer region table output of ROSE_Main and map genes to it
#will create two outputs a gene mapped region table where each row is an enhancer
#and a gene table where each row is a gene
#does this by default for super-enhancers only
import ROSE_utils


def getDistance(locus,othreLocus):
    distance = (locus.start()-othreLocus.start() + locus.end()-othreLocus.end())/2
    return abs(distance)

#==================================================================
#====================MAPPING GENES TO DHSS====================
#==================================================================
def mapEnhancerToGene(annotFile,enhancerFile,transcribedFile='',uniqueGenes=True,searchWindow =100000):
    startDict = ROSE_utils.makeStartDict(annotFile)
    enhancerTable = ROSE_utils.parseTable(enhancerFile,'\t')
    if len(transcribedFile) > 0:
        transcribedTable = ROSE_utils.parseTable(transcribedFile,'\t')
        transcribedGenes = [line[1] for line in transcribedTable]
    else:
        transcribedGenes = startDict.keys()
    print('MAKING TSS COLLECTION')
    tssLoci = []
    for geneID in transcribedGenes:
        tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,0,0))
    #this turns the tssLoci list into a LocusCollection
    #50 is the internal parameter for LocusCollection and doesn't really matter
    tssCollection = ROSE_utils.LocusCollection(tssLoci,50)
    enhancerToGeneTable = [enhancerTable[0]+['GENES','DISTANCE']]

    for line in tqdm(enhancerTable):
        if line[0][0] =='#' or line[0][0] == 'R':
            continue
        enhancerLocus = ROSE_utils.Locus(line[1],line[2],line[3],'.',line[0])

        proximalLoci = tssCollection.getOverlap(ROSE_utils.makeSearchLocus(enhancerLocus,searchWindow,searchWindow),'both')           
        proximalGenes = []
        for proxLocus in proximalLoci:
            proximalGenes.append([proxLocus.ID(),getDistance(enhancerLocus,proxLocus)])

        #NOW WRITE THE ROW FOR THE DHS TABLE
        newEnhancerLine = list(line)
        newEnhancerLines = [newEnhancerLine + genes for genes in proximalGenes]

        enhancerToGeneTable.extend(newEnhancerLines)

    return enhancerToGeneTable

#==================================================================
#=========================MAIN METHOD==============================
#==================================================================

def main():
    from optparse import OptionParser
    usage = "usage: %prog [options] -g [GENOME] -i [INPUT_DHS_FILE]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-i","--i", dest="input",nargs = 1, default=None,
                      help = "Enter a ROSE ranked enhancer or super-enhancer file")
    parser.add_option("-g","--genome", dest="genome",nargs = 1, default=None,
                      help = "Enter the genome build (MM10,MM9,MM8,HG18,HG19,HG38)")

    #optional flags
    parser.add_option("-l","--list", dest="geneList",nargs = 1, default=None,
                      help = "Enter a gene list to filter through")
    parser.add_option("-o","--out", dest="out",nargs = 1, default=None,
                      help = "Enter an output folder. Default will be same folder as input file")
    parser.add_option("-w","--window", dest="window",nargs = 1, default=100000,
                      help = "Enter a search distance for genes. Default is 100,000bp")
    #RETRIEVING FLAGS
    (options,args) = parser.parse_args()


    if not options.input or not options.genome:
        parser.print_help()
        exit()

    #GETTING THE INPUT
    enhancerFile = options.input
    window = int(options.window)

    #making the out folder if it doesn't exist
    if options.out:
        outFolder = ROSE_utils.formatFolder(options.out,True)
    else:
        outFolder = '/'.join(enhancerFile.split('/')[0:-1]) + '/'


    #GETTING THE GENOME
    genome = options.genome
    print('USING %s AS THE GENOME' % genome)

    #GETTING THE CORRECT ANNOT FILE
    cwd = os.getcwd()
    genomeDict = {
        'HG18':'%s/annotation/hg18_refseq.ucsc' % (cwd),
        'MM9': '%s/annotation/mm9_refseq.ucsc' % (cwd),
        'HG19':'%s/annotation/hg19_refseq.ucsc' % (cwd),
	    'HG38':'%s/annotation/hg38_refseq.ucsc' % (cwd),
        'MM8': '%s/annotation/mm8_refseq.ucsc' % (cwd),
        'MM10':'%s/annotation/mm10_refseq.ucsc' % (cwd),
        }

    annotFile = genomeDict[genome.upper()]

    #GETTING THE TRANSCRIBED LIST
    if options.geneList:
        transcribedFile = options.geneList
    else:
        transcribedFile = ''

    enhancerToGeneTable = mapEnhancerToGene(annotFile,enhancerFile,transcribedFile,True,window )

    #Writing enhancer output
    enhancerFileName = enhancerFile.split('/')[-1].split('.')[0]

    if window != 100000:
        #writing the enhancer table
        out1 = '%s%s_DHS_TO_GENE_%sKB.txt' % (outFolder,enhancerFileName,window/1000)
        ROSE_utils.unParseTable(enhancerToGeneTable,out1,'\t')
    else:
        #writing the enhancer table
        out1 = '%s%s_DHS_TO_GENE.txt' % (outFolder,enhancerFileName)
        ROSE_utils.unParseTable(enhancerToGeneTable,out1,'\t')
if __name__ == "__main__":
    main()
