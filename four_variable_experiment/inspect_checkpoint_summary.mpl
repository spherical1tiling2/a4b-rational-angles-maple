restart:
read "three_to_two_fullgrid_small.m":
printf("SMALL_PROCESSED=%a\n",processedCases):
printf("SMALL_RANK_DEFICIENT=%a\n",rankDeficientCases):
printf("SMALL_DEFERRED=%a\n",deferredCases):
printf("SMALL_CASE_RESULT_COUNT=%d\n",numelems([indices(caseResults)])):
quit:
