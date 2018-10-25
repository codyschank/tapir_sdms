# tapir_sdms

**Order of steps to complete modeling:**
1. Run behindTheScenes_predictorPrep.R. This will do the initiall legwork to download and process environmental layers.
2. Run predictorPrep.r on the results from step 1. This will do the final preparation of environmental layers before modelling. predictorPrep_4km.R is a specialized case of this prep.
3. Run presencePrep.R on your presence data.
4. Run either chapter 3 or chapter 4 model loop.

**Next steps:**
1. Make changes to the code so that scripts are more easily applied to other species and regions
