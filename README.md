# newCET_uai2023
Matlab code for UAI2023 paper "Establishing Markov equivalence in cyclic directed graphs" (Claassen & Mooij, UAI2023) PMLR 216:433-442
NOTE: The 'cg_to_cpag_new' algorithm, stage (f) includes the adjustment to Algorithm 2 resulting from the correction to rule (iv) in Theorem 1, see new article version on arXiv for details.

Main scripts:
- cg_to_cpag_new.m 	= convert (cyclic) directed graph via CMAG into CPAG representation (combines Algorithms 1+2 from the main article)
- cpag_from_cg_org.m 	= original version based on d-separation test from 'Discovering cyclic causal structure' (Richardson, 1996)
- mk_random_cg.m 	= create random cyclic graph over N nodes (with configurable density/cycles)
- csep.m 		= d-separation test in cyclic graph
- get_scc_an.m		= partition cyclic graph into strongly connected components + ancestral matrix

Experimental evaluation:
- UAI_run_CPAG_test.m 	= run batch test for random cyclic graph to CPAG algorithms (org+new)
- analRes_UAIrun.m 	= convert results batch run into figures 4 (main) + 5 (supplement)

Tools:
- myunion/intersect/setdiff = faster versions of built in set manipulation routines
- draw_cpmag.m 		= procedure to visualise (cyclic) directed graph
