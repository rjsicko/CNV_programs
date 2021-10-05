example internal cumulative database. It is labeled "cleaned" for the following reasons:

- we only included analyses clustered on our samples, so there is less noise in the data.
- HapMap samples were excluded (they generally performed different than our samples on the array; therefore, the CNV calls we generate with them might be an artifact of us clustering on our data)
- any sample with a LRR deviation > 0.2 (overall measure of how messy a sample is) excluded

