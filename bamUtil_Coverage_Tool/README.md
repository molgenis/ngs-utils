This tool is calculating coverage per position of pathogenic and likely pathogenic varaints from the Managed Variant List (MVL) (UMCG-MVL_VKGLconsensusMVL-2017-11-27).

run MVL_Coverage_Tool.sh to execute

required input is a folder containing BAM files.
    -s|                 locaion of folder with bam files
    -t|			Give loaction of tmp folder for inbetween and final files
optional input
	-b|--bamstat		path to the bam stats (default: $EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/bamUtil-1.0.14/bin/bam stats)
	-m|--mvl	        location of MVL_per_bp_unique.txt (default:$EBROOTNGSMINUTILS/bamUtil_Coverage_Tool//MVL_per_bp_unique.txt)
	-x|--mvltabix       location of MVL_per_bp_unique_tabix.txt (default: $EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/MVL_per_bp_unique_tabix.txt)
	-q|--qcmvldir       location of MerdgeQCofMVL_file.py necessary for final txt generation (default:$EBROOTNGSMINUTILS/bamUtil_Coverage_Tool/)

 -b     The bam stats tool is located here $EBROOTNGSMINUTILS/bamUtil_Coverage_Tool//bamUtil-1.0.14/bin/bam stats.

-mx     you can make your own position list or MVL, make sure you also make a tabix.txt of this file. even if you give regions, the output will be per bp.

-t      location for in-between and final files.

-q      a python scrip generating the 4 txt files.
        Avrg_SD_moreThan30x20x10x.txt: all postions with average coverage, SD, percenateg more than 30x, 20x and 10x.
        Avrg_SD_moretha20x10x_Molgenis.txt: like Avrg_SD_moreThan30x20x10x.txt but fit for Molgenis96
        MVL_outliers.txt: all the outliers of the input sample set, if coverage is above or below 3x the SD from the average, it is in this txt file.
        PosBelow_30x_10x_1x.txt: all position with an average coverage below 30x, 10x and 1x are in this txt file.
	
