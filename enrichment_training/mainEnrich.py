import os

#= perl compute_enrichment.pl [namecode.txt] [outputfolder] 1
#= 1: upper cumulative
#= 2: lower comulative
#cmd = 'perl compute_enrichment.pl ../output/preXLOC_GO/namecode.txt ../output/preXLOC_GO/output_comput_enrich 1';os.system(cmd)
cmd = 'perl compute_enrichment.pl output/top5scoredTG_pathways/namecode.txt output/top5scoredTG_pathways/output_comput_enrich 1';os.system(cmd)
